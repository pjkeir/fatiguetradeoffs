# Functions to calcluate joint positions, segment COM, and whole body COM location in a specific axis direction

import math

def jointpositions(segment_angles, segment_lengths):

    # Horizontal location of the joints (x-position); right endpoint is (0, 0)
    x1 = segment_lengths[0] * math.cos(math.radians(segment_angles[0]))                # Right wrist
    x2 = x1 + segment_lengths[1] * math.cos(math.radians(segment_angles[1]))           # Right elbow
    x3 = x2 + segment_lengths[2] * math.cos(math.radians(segment_angles[2]))           # Right shoulder
    x4 = x3 + segment_lengths[3] * math.cos(math.radians(segment_angles[6]))           # Low back
    x5 = x4 + segment_lengths[4] * math.cos(math.radians(segment_angles[7]))           # Knee
    x6 = x5 + segment_lengths[5] * math.cos(math.radians(segment_angles[8]))           # Ankle
    x30 = x3 + segment_lengths[2] * math.cos(math.radians(segment_angles[5])) * -1     # Left elbow
    x20 = x30 + segment_lengths[1] * math.cos(math.radians(segment_angles[4])) * -1    # Left wrist
    x10 = x20 + segment_lengths[0] * math.cos(math.radians(segment_angles[3])) * -1    # Left endpoint
    x_val = [x1, x2, x3, x4, x5, x6, x30, x20, x10]

    # Vertical location of the joints (y-position); right endpoint is (0, 0)
    y1 = segment_lengths[0] * math.sin(math.radians(segment_angles[0])) * -1           # Right wrist
    y2 = y1 + segment_lengths[1] * math.sin(math.radians(segment_angles[1])) * -1      # Right elbow
    y3 = y2 + segment_lengths[2] * math.sin(math.radians(segment_angles[2])) * -1      # Right shoulder
    y4 = y3 + segment_lengths[3] * math.sin(math.radians(segment_angles[6])) * -1      # Low back
    y5 = y4 + segment_lengths[4] * math.sin(math.radians(segment_angles[7])) * -1      # Knee
    y6 = y5 + segment_lengths[5] * math.sin(math.radians(segment_angles[8])) * -1      # Ankle
    y30 = y3 + segment_lengths[2] * math.sin(math.radians(segment_angles[5]))          # Left elbow
    y20 = y30 + segment_lengths[1] * math.sin(math.radians(segment_angles[4]))         # Left wrist
    y10 = y20 + segment_lengths[0] * math.sin(math.radians(segment_angles[3]))         # Left endpoint
    y_val = [y1, y2, y3, y4, y5, y6, y30, y20, y10]

    positions = [x_val, y_val]

    return positions

def segmentcenterofmass(positions, segment_com_loc):

    hand_com = positions[0] * ((100 - segment_com_loc[0]) / 100)
    forearm_com = positions[0] + (positions[1] - positions[0]) * ((100 - segment_com_loc[1]) / 100)
    upperarm_com = positions[1] + (positions[2] - positions[1]) * ((100 - segment_com_loc[2]) / 100)
    torso_com = positions[2] + (positions[3] - positions[2]) * ((100 - segment_com_loc[3]) / 100)
    thigh_com = positions[3] + (positions[4] - positions[3]) * ((100 - segment_com_loc[4]) / 100)
    shank_com = positions[4] + (positions[5] - positions[4]) * ((100 - segment_com_loc[5]) / 100)
    left_upperarm_com = positions[6] + (positions[2] - positions[6]) * ((100 - segment_com_loc[2]) / 100)
    left_forearm_com = positions[7] + (positions[6] - positions[7]) * ((100 - segment_com_loc[1]) / 100)
    left_hand_com = positions[8] + (positions[7] - positions[8]) * ((100 - segment_com_loc[0]) / 100)

    return[hand_com, forearm_com, upperarm_com, torso_com, thigh_com, shank_com, left_upperarm_com, left_forearm_com, left_hand_com]

def bodycenterofmass(positions, segment_mass, segment_com_loc):

    segment_com = segmentcenterofmass(positions, segment_com_loc)

    summedmass = (segment_mass[0] * 2       # Hands (2)
                  + segment_mass[1] * 2     # Forearms (2)
                  + segment_mass[2] * 2     # Upperarm (2)
                  + segment_mass[3]         # Trunk (lumped)
                  + segment_mass[4]         # Thigh (lumped)
                  + segment_mass[5])        # Shank (lumped)

    body_com = (segment_com[0] * (segment_mass[0] / summedmass)
                  + segment_com[1] * (segment_mass[1] / summedmass)
                  + segment_com[2] * (segment_mass[2] / summedmass)
                  + segment_com[3] * (segment_mass[3] / summedmass)
                  + segment_com[4] * (segment_mass[4] / summedmass)
                  + segment_com[5] * (segment_mass[5] / summedmass)
                  + segment_com[6] * (segment_mass[2] / summedmass)
                  + segment_com[7] * (segment_mass[1] / summedmass)
                  + segment_com[8] * (segment_mass[0] / summedmass))

    return body_com
