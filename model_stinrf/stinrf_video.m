function out = stinrf_video(M, framerate, name, varargin)
% STINRF_VIDEO  Write one stimulus as a video, from the matrix a figure uses.
%
%   The export scripts call this for every stimulus their figure shows, in
%   the same run and from the same matrix as the panel, so the video and the
%   figure cannot end up showing two different stimuli.
%
%   A 1D stimulus (rows are frames, columns are space) is drawn as vertical
%   stripes filling a square frame, the way it looks on a screen; a 2D one
%   (rows x columns x frames) is written as it is. Both are enlarged by a
%   whole number, pixel for pixel, to at least MinWidth across.
%
%   Grey levels are luminance, not stretched to fill the range: Clim(1) is
%   black and Clim(2) white. With Clim = [0 2*meanL] the mean luminance is
%   mid grey and a low-contrast stimulus looks low-contrast. Values outside
%   Clim are clipped, with a warning.
%
%   Playback is real time. A screen shows at most about 120 frames a second,
%   so a stimulus sampled faster is played one frame in k, k being the
%   smallest whole number that brings the rate down to MaxFps: 120 fps plays
%   every frame, 200 fps plays at 100 and 1000 fps at 111.1. A Speed below 1
%   plays it in slow motion instead; at Speed = 0.12 a 1000 fps stimulus
%   shows every one of its frames. A 120 fps stimulus that changes on every
%   frame, such as reverse phi, needs a 120 Hz screen to be seen as built.
%
%   Where MATLAB can write MPEG-4 (Windows and macOS) it does. On Linux it
%   cannot, so the frames go through ffmpeg; with neither, the video is
%   written as a Motion JPEG AVI.
%
% Inputs
%   M (double) T x N (1D) or H x W x T (2D) luminance
%   framerate (double) frame rate the stimulus was built at [fps]
%   name (char) file name without extension, e.g. 'SecondOrder_1_Leftward';
%       each run of characters other than letters, digits, '-' and '.'
%       becomes one '_'
%
% Name-value parameters:
%   'Clim' ([lo hi]) luminance shown as black and white. Default: the
%       stimulus's own minimum and maximum
%   'Speed' (double) playback speed, 1 = real time. Default 1
%   'MaxFps' (double) highest frame rate written. Default 120
%   'MinWidth' (int) smallest width of the video [px]. Default 600
%   'OutDir' (char) folder. Defaults to the STINRF_VIDEOS environment
%       variable if set, otherwise a videos folder beside this file.
%       STINRF_VIDEOS=off writes nothing
%
% Outputs
%   out (char) path of the file written, '' when videos are off

    p = inputParser;
    p.addParameter('Clim', []);
    p.addParameter('Speed', 1);
    p.addParameter('MaxFps', 120);
    p.addParameter('MinWidth', 600);
    p.addParameter('OutDir', default_video_dir());
    p.parse(varargin{:});
    o = p.Results;

    out = '';
    if strcmpi(o.OutDir, 'off'), return; end
    if ~exist(o.OutDir, 'dir'), mkdir(o.OutDir); end
    name = regexprep(char(name), '[^A-Za-z0-9.\-]+', '_');

    is1D = ismatrix(M);
    if is1D, T = size(M, 1); W = size(M, 2); H = W;
    else,    T = size(M, 3); W = size(M, 2); H = size(M, 1);
    end

    % Whole-number decimation keeps the frames evenly spaced in time;
    % resampling to exactly 120 fps would repeat some frames and not others.
    k   = max(1, ceil(framerate * o.Speed / o.MaxFps - 1e-9));
    idx = 1:k:T;
    fps = framerate * o.Speed / k;

    if isempty(o.Clim), o.Clim = [min(M(:)), max(M(:))]; end
    lo = o.Clim(1);  hi = o.Clim(2);
    tol = 1e-9 * max(1, abs(hi - lo));
    clipped = mean(M(:) < lo - tol | M(:) > hi + tol);
    if clipped > 0
        warning('stinrf_video:clip', ['%s: %.2f%% of the values fall ' ...
                'outside [%g %g] and are clipped'], name, 100*clipped, lo, hi);
    end

    s  = max(1, ceil(o.MinWidth / W));
    Hv = H*s + mod(H*s, 2);            % H.264 wants even dimensions
    Wv = W*s + mod(W*s, 2);
    frame = @(i) grey(M, idx(i), is1D, H, lo, hi, s, Hv, Wv);

    profiles = {VideoWriter.getProfiles().Name};
    if any(strcmp(profiles, 'MPEG-4'))
        out = fullfile(o.OutDir, [name '.mp4']);
        v = VideoWriter(out, 'MPEG-4');
        v.FrameRate = fps;  v.Quality = 100;
        open(v);
        for i = 1:numel(idx), writeVideo(v, repmat(frame(i), 1, 1, 3)); end
        close(v);
    elseif has_ffmpeg()
        out = fullfile(o.OutDir, [name '.mp4']);
        raw = [tempname '.gray'];
        fid = fopen(raw, 'w');
        % ffmpeg reads raw frames row by row; MATLAB stores them column by
        % column, hence the transpose.
        for i = 1:numel(idx), fwrite(fid, frame(i).', 'uint8'); end
        fclose(fid);
        [n, d] = rat(fps, 1e-9);
        [st, msg] = system(sprintf([ffmpeg_prefix() 'ffmpeg -nostdin -y -loglevel error ' ...
            '-f rawvideo -pix_fmt gray -s %dx%d -r %d/%d -i "%s" ' ...
            '-c:v libx264 -preset slow -crf 10 -pix_fmt yuv420p "%s"'], ...
            Wv, Hv, n, d, raw, out));
        delete(raw);
        if st ~= 0, error('stinrf_video: ffmpeg failed on %s:\n%s', name, msg); end
    else
        out = fullfile(o.OutDir, [name '.avi']);
        warning('stinrf_video:avi', ['%s: this MATLAB cannot write MPEG-4 ' ...
                'and ffmpeg was not found, so it is a Motion JPEG AVI'], name);
        v = VideoWriter(out, 'Motion JPEG AVI');
        v.FrameRate = fps;  v.Quality = 100;
        open(v);
        for i = 1:numel(idx), writeVideo(v, frame(i)); end
        close(v);
    end
    fprintf('  video %s  (%d frames at %.4g fps, %dx%d)\n', ...
            out, numel(idx), fps, Wv, Hv);
end


function G = grey(M, t, is1D, H, lo, hi, s, Hv, Wv)
% One frame as 8-bit grey, enlarged s times and padded to Hv x Wv.
    if is1D, G = repmat(M(t, :), H, 1);
    else,    G = M(:, :, t);
    end
    G = uint8(255 * min(max((G - lo) / (hi - lo), 0), 1));
    G = repelem(G, s, s);
    G(end+1:Hv, :) = repmat(G(end, :), Hv - size(G, 1), 1);
    G(:, end+1:Wv) = repmat(G(:, end), 1, Wv - size(G, 2));
end


function tf = has_ffmpeg()
    [st, ~] = system([ffmpeg_prefix() 'ffmpeg -version']);
    tf = st == 0;
end


function p = ffmpeg_prefix()
% MATLAB on Linux puts its own libraries first on LD_LIBRARY_PATH, and a
% system ffmpeg started from inside it can then fail to load. Clearing the
% variable for that one command lets it find its own.
    if isunix && ~ismac, p = 'LD_LIBRARY_PATH= '; else, p = ''; end
end


function d = default_video_dir()
% Where the videos go unless told otherwise; mirrors default_export_dir.
    d = getenv('STINRF_VIDEOS');
    if isempty(d)
        d = fullfile(inme_output(), 'videos');
    end
end
