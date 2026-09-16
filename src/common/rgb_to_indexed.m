function [ind, cmap] = rgb_to_indexed(rgb, n_colors)
%RGB_TO_INDEXED  Quantise a truecolour frame to an indexed image, portably.
%
%   [ind, cmap] = RGB_TO_INDEXED(rgb)
%   [ind, cmap] = RGB_TO_INDEXED(rgb, n_colors)
%
%   MATLAB accepts rgb2ind(RGB, N) to request an N-entry colourmap, which is
%   what keeps an animated GIF small. GNU Octave implements only the
%   single-argument form, and that returns one colourmap entry per distinct
%   colour in the frame. On an anti-aliased plot that is tens of thousands of
%   entries per frame and produces GIFs two orders of magnitude too large.
%
%   When the explicit form is unavailable this falls back to snapping each
%   channel onto a uniform cube first, bounding the palette before conversion.

    if nargin < 2 || isempty(n_colors)
        n_colors = 256;
    end

    try
        [ind, cmap] = rgb2ind(rgb, n_colors);
    catch
        % levels^3 <= n_colors, so 6 levels gives a 216-colour cube.
        levels = max(2, floor(n_colors^(1/3)));
        x = double(rgb);
        if ~isinteger(rgb) && max(x(:)) <= 1
            x = x * 255;
        end
        q = round(x / 255 * (levels-1)) / (levels-1);
        [ind, cmap] = rgb2ind(q);
    end
end
