% clamps the given coordinate to the range passed.
% coord: the given pair of coordinates
% range: [TOP RIGHT BOTTOM LEFT]
function clamped = clamp_to_range(coord, range)
    clamped(1) = max(range(4), min(range(2), coord(1)));
    clamped(2) = max(range(1), min(range(3), coord(2)));
end