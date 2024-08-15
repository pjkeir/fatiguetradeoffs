# Inverse Dynamics to output internal joint moments given posture, model, and external load
# Two slightly different functions, both give equivalent results but #1 is 2x faster
# inversedynamics1: uses trigonometric equations
# inversedynamics2: uses cross-products


import math
import numpy as np
import posture

def inversedynamics1(theta, segment_mass, segment_lengths, segment_com_loc, load, loadangle):

    # Convert segment lengths from cm to m
    segment_lengths = np.divide(segment_lengths, 100)

    # Right wrist (assume load acts at COM)
    t0 = math.radians(theta[0])
    fx0 = 9.81 * load * math.cos(loadangle)
    fy0 = -9.81 * segment_mass[0] + load * 9.81 * math.sin(loadangle)
    moment0 = fx0 * segment_lengths[0] * (segment_com_loc[0] / 100) * math.sin(t0) \
              + fy0 * segment_lengths[0] * (segment_com_loc[0] / 100) * math.cos(t0)

    # Right elbow
    t1 = math.radians(theta[1])
    fx1 = fx0
    fy1 = fy0 + -9.81 * segment_mass[1]
    moment1 = moment0 \
              + fx0 * segment_lengths[1] * math.sin(t1) \
              + fy0 * segment_lengths[1] * math.cos(t1) \
              + -9.81 * segment_mass[1] * segment_lengths[1] * (segment_com_loc[1] / 100) * math.cos(t1)

    # Right shoulder
    t2 = math.radians(theta[2])
    fx2 = fx1
    fy2 = fy1 + -9.81 * segment_mass[2]
    moment2 = moment1 \
              + fx1 * segment_lengths[2] * math.sin(t2) \
              + fy1 * segment_lengths[2] * math.cos(t2) \
              + -9.81 * segment_mass[2] * segment_lengths[2] * (segment_com_loc[2] / 100) * math.cos(t2)

    # Left wrist
    t3 = math.radians(theta[3])
    fx3 = 0
    fy3 = -9.81 * segment_mass[0]
    moment3 = fx3 * segment_lengths[0] * (segment_com_loc[0] / 100) * math.sin(t3) \
              + fy3 * segment_lengths[0] * (segment_com_loc[0] / 100) * math.cos(t3)

    # Left elbow
    t4 = math.radians(theta[4])
    fx4 = fx3
    fy4 = fy3 + -9.81 * segment_mass[1]
    moment4 = moment3 \
              + fx3 * segment_lengths[1] * math.sin(t4) \
              + fy3 * segment_lengths[1] * math.cos(t4) \
              + -9.81 * segment_mass[1] * segment_lengths[1] * (segment_com_loc[1] / 100) * math.cos(t4)

    # Left shoulder
    t5 = math.radians(theta[5])
    fx5 = fx4
    fy5 = fy4 + -9.81 * segment_mass[2]
    moment5 = moment4 \
              + fx4 * segment_lengths[2] * math.sin(t5) \
              + fy4 * segment_lengths[2] * math.cos(t5) \
              + -9.81 * segment_mass[2] * segment_lengths[2] * (segment_com_loc[2] / 100) * math.cos(t5)

    # Low back
    t6 = math.radians(theta[6])
    fx6 = fx2 + fx5
    fy6 = fy2 + fy5 + -9.81 * segment_mass[3]
    moment6 = moment2 + moment5 \
              + (fx2 + fx5) * segment_lengths[3] * math.sin(t6) \
              + (fy2 + fy5) * segment_lengths[3] * math.cos(t6) \
              + -9.81 * segment_mass[3] * segment_lengths[3] * (segment_com_loc[3] / 100) * math.cos(t6)

    # Knee
    t7 = math.radians(theta[7])
    fx7 = fx6
    fy7 = fy6 + -9.81 * segment_mass[4]
    moment7 = moment6 \
              + fx6 * segment_lengths[4] * math.sin(t7) \
              + fy6 * segment_lengths[4] * math.cos(t7) \
              + -9.81 * segment_mass[4] * segment_lengths[4] * (segment_com_loc[4] / 100) * math.cos(t7)

    # Ankle
    t8 = math.radians(theta[8])
    fx8 = fx7
    fy8 = fy7 + -9.81 * segment_mass[5]
    moment8 = moment7 \
              + fx7 * segment_lengths[5] * math.sin(t8) \
              + fy7 * segment_lengths[5] * math.cos(t8) \
              + -9.81 * segment_mass[5] * segment_lengths[5] * (segment_com_loc[5] / 100) * math.cos(t8)

    return [ moment0, moment1, moment2, moment3, moment4, moment5, moment6, moment7, moment8 ]

def inversedynamics2(theta, segment_mass, segment_lengths, segment_com_loc, load, loadangle):

    # Convert segment lengths from cm to m
    segment_lengths = np.divide(segment_lengths, 100)

    # Determining joint positions, COMx, and COMy locations
    positions = posture.jointpositions(theta, segment_lengths)
    segmentcom_x = posture.segmentcenterofmass(positions[0], segment_com_loc)
    segmentcom_y = posture.segmentcenterofmass(positions[1], segment_com_loc)

    # Right wrist(assume load acts at COM)
    t0 = math.radians(theta[0])
    handload = 9.81 * load * np.array([math.cos(loadangle), math.sin(loadangle)])
    hand_wt = 9.81 * segment_mass[0] * np.array([0, -1])
    handcom_pos_wrist = np.array([segmentcom_x[0], segmentcom_y[0]]) - np.array([positions[0][0], positions[1][0]])
    internalmoment0 = -1 * (np.cross(handcom_pos_wrist, handload) + np.cross(handcom_pos_wrist, hand_wt))
    jrf0 = -1 * (handload + hand_wt)

    # Right elbow
    t1 = math.radians(theta[1])
    forearm_wt = 9.81 * segment_mass[1] * np.array([0, -1])
    forearmcom_pos_elbow = np.array([segmentcom_x[1], segmentcom_y[1]]) - np.array([positions[0][1], positions[1][1]])
    wrist_pos_elbow = np.array([positions[0][0], positions[1][0]]) - np.array([positions[0][1], positions[1][1]])
    internalmoment1 = -1 * (-1 * internalmoment0 + np.cross(wrist_pos_elbow, -1 * jrf0) + np.cross(forearmcom_pos_elbow, forearm_wt))
    jrf1 = -1 * (forearm_wt + -1 * jrf0)

    # Right shoulder
    t2 = math.radians(theta[2])
    upperarm_wt = 9.81 * segment_mass[2] * np.array([0, -1])
    upperarmcom_pos_shoulder = np.array([segmentcom_x[2], segmentcom_y[2]]) - np.array([positions[0][2], positions[1][2]])
    elbow_pos_shoulder =  np.array([positions[0][1], positions[1][1]]) - np.array([positions[0][2], positions[1][2]])
    internalmoment2 = -1 * (-1 * internalmoment1 + np.cross(elbow_pos_shoulder, -1 * jrf1) + np.cross(upperarmcom_pos_shoulder, upperarm_wt))
    jrf2 = -1 * (upperarm_wt + -1 * jrf1)

    # Left wrist
    t3 = math.radians(theta[3])
    lefthand_wt = 9.81 * segment_mass[0] * np.array([0, -1])
    lefthandcom_pos_wrist = np.array([segmentcom_x[8], segmentcom_y[8]]) - np.array([positions[0][7], positions[1][7]])
    internalmoment3 = -1 * (np.cross(lefthandcom_pos_wrist, lefthand_wt))
    jrf3 = -1 * (lefthand_wt)

    # Left elbow
    t4 = math.radians(theta[4])
    leftforearm_wt = 9.81 * segment_mass[1] * np.array([0, -1])
    leftforearmcom_pos_elbow = np.array([segmentcom_x[7], segmentcom_y[7]]) - np.array([positions[0][6], positions[1][6]])
    leftwrist_pos_elbow = np.array([positions[0][7], positions[1][7]]) - np.array([positions[0][6], positions[1][6]])
    internalmoment4 = -1 * (-1 * internalmoment3 + np.cross(leftwrist_pos_elbow, -1 * jrf3) + np.cross(leftforearmcom_pos_elbow, leftforearm_wt))
    jrf4 = -1 * (leftforearm_wt + -1 * jrf3)

    # Left shoulder
    t5 = math.radians(theta[5])
    leftupperarm_wt = 9.81 * segment_mass[2] * np.array([0, -1])
    leftupperarmcom_pos_shoulder = np.array([segmentcom_x[6], segmentcom_y[6]]) - np.array([positions[0][2], positions[1][2]])
    leftelbow_pos_shoulder = np.array([positions[0][6], positions[1][6]]) - np.array([positions[0][2], positions[1][2]])
    internalmoment5 = -1 * (-1 * internalmoment4 + np.cross(leftelbow_pos_shoulder, -1 * jrf4) + np.cross(leftupperarmcom_pos_shoulder, leftupperarm_wt))
    jrf5 = -1 * (leftupperarm_wt + -1 * jrf4)

    # Low back
    t6 = math.radians(theta[6])
    back_wt = 9.81 * segment_mass[3] * np.array([0, -1])
    backcom_pos_backjoint = np.array([segmentcom_x[3], segmentcom_y[3]]) - np.array([positions[0][3], positions[1][3]])
    shoulder_pos_backjoint = np.array([positions[0][2], positions[1][2]]) - np.array([positions[0][3], positions[1][3]])
    internalmoment6 = -1 * (-1 * (internalmoment2 + internalmoment5) + np.cross(shoulder_pos_backjoint, -1 * (jrf2 + jrf5)) + np.cross(backcom_pos_backjoint, back_wt))
    jrf6 = -1 * (back_wt + -1 * (jrf2 + jrf5))

    # Knee
    t7 = math.radians(theta[7])
    thigh_wt = 9.81 * segment_mass[4] * np.array([0, -1])
    thighcom_pos_knee = np.array([segmentcom_x[4], segmentcom_y[4]]) - np.array([positions[0][4], positions[1][4]])
    backjoint_pos_knee = np.array([positions[0][3], positions[1][3]]) - np.array([positions[0][4], positions[1][4]])
    internalmoment7 = -1 * (-1 * internalmoment6 + np.cross(backjoint_pos_knee, -1 * jrf6) + np.cross(thighcom_pos_knee, thigh_wt))
    jrf7 = -1 * (thigh_wt + -1 * jrf6)

    # Ankle
    t8 = math.radians(theta[8])
    shank_wt = 9.81 * segment_mass[5] * np.array([0, -1])
    shankcom_pos_ankle = np.array([segmentcom_x[5], segmentcom_y[5]]) - np.array([positions[0][5], positions[1][5]])
    knee_pos_ankle = np.array([positions[0][4], positions[1][4]]) - np.array([positions[0][5], positions[1][5]])
    internalmoment8 = -1 * (-1 * internalmoment7 + np.cross(knee_pos_ankle, -1 * jrf7) + np.cross(shankcom_pos_ankle, shank_wt))
    jrf8 = -1 * (shank_wt + -1 * jrf7)

    return [internalmoment0, internalmoment1, internalmoment2, internalmoment3, internalmoment4,
            internalmoment5, internalmoment6, internalmoment7, internalmoment8]
