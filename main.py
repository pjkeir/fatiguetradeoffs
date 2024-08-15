"""
MAIN SCRIPT:

Given input simulation parameters for a one-arm pulling task, this script will randomly sample
anthropometrics, calculate segment parameters, predict posture based on minimizing joint moments, and
calculate joint moments when that posture is perturbed across the upper extremity joints to
simulate "subtle" fatigue-induced changes in joint kinematics.

Created by: Kevin Kos
Edited by: Daanish Mulla
Project Supervisor: Dr. Peter Keir

"""

# Libraries
import datetime
import math
import numpy as np
import pandas as pd
from scipy.optimize import minimize
import sys


# Modules
import moments
import posture
import segments
import visualize


"""
Simulation Parameters
"""
np.random.seed(2022)                               # Setting  seed for reproducing results
nSubjects = 500                                    # Number of subjects simulated
joints_perturbed = ['Wrist', 'Elbow', 'Shoulder']  # Joints perturbed
perturbation = 5                                   # Magnitude of +/- joint angle perturbation (degrees)
reach_heights_simulated = [80, 100, 120, 140]      # Reach heights
sexs = ['M', 'F']                                  # Sex
load = 50/9.81                                     # Load magnitude (default set as 50 N)
loadangle = math.radians(180)                      # Flat, pulling (in the x direction); 0 = push; 180 = pull
constrainFoot = False                              # Constrain horizontal ankle location in absolute distance
constrainDistance = 50                             # Horizontal ankle location fixed (in cm)
constrainRelative = True                           # Constrain horizontal ankle location relative to arm length
constrainRelativePercent = 60                      # Relative percentage of arm length
failedsolutions = 0                                # Counter for number of failed simulations

"""
Visualization Parameters
"""
visualizeInitialguess = False                      # Visualizes the initial posture guess
visualizeInitialsolution = False                   # Visualizes the initial optimal solution for a given virtual subject / task
visualizechange = False                            # Visualizes the change in solution with joint perturbations

"""
Check
"""
evalaverage = False                                # An indicator on whether to only look at the average individual (for re-creating Figure 1)
if evalaverage is True:
    nSubjects = 1
    visualizeInitialsolution = True


"""
Placeholder Array for Data Storage
"""
simulationresults = []                             # Storing results of simulation for statistical analysis

# Storing data for checking convergence (Males and Females separate)
optimal_rwmoments_M = np.zeros((nSubjects, len(reach_heights_simulated)))
optimal_remoments_M = np.zeros((nSubjects, len(reach_heights_simulated)))
optimal_rsmoments_M = np.zeros((nSubjects, len(reach_heights_simulated)))
optimal_bmoments_M = np.zeros((nSubjects, len(reach_heights_simulated)))
optimal_rwmoments_M_cummean = np.zeros((nSubjects, len(reach_heights_simulated)))
optimal_remoments_M_cummean = np.zeros((nSubjects, len(reach_heights_simulated)))
optimal_rsmoments_M_cummean = np.zeros((nSubjects, len(reach_heights_simulated)))
optimal_bmoments_M_cummean = np.zeros((nSubjects, len(reach_heights_simulated)))
optimal_rwmoments_M_cumstd = np.zeros((nSubjects, len(reach_heights_simulated)))
optimal_remoments_M_cumstd = np.zeros((nSubjects, len(reach_heights_simulated)))
optimal_rsmoments_M_cumstd = np.zeros((nSubjects, len(reach_heights_simulated)))
optimal_bmoments_M_cumstd = np.zeros((nSubjects, len(reach_heights_simulated)))
optimal_rwmoments_F = np.zeros((nSubjects, len(reach_heights_simulated)))
optimal_remoments_F = np.zeros((nSubjects, len(reach_heights_simulated)))
optimal_rsmoments_F = np.zeros((nSubjects, len(reach_heights_simulated)))
optimal_bmoments_F = np.zeros((nSubjects, len(reach_heights_simulated)))
optimal_rwmoments_F_cummean = np.zeros((nSubjects, len(reach_heights_simulated)))
optimal_remoments_F_cummean = np.zeros((nSubjects, len(reach_heights_simulated)))
optimal_rsmoments_F_cummean = np.zeros((nSubjects, len(reach_heights_simulated)))
optimal_bmoments_F_cummean = np.zeros((nSubjects, len(reach_heights_simulated)))
optimal_rwmoments_F_cumstd = np.zeros((nSubjects, len(reach_heights_simulated)))
optimal_remoments_F_cumstd = np.zeros((nSubjects, len(reach_heights_simulated)))
optimal_rsmoments_F_cumstd = np.zeros((nSubjects, len(reach_heights_simulated)))
optimal_bmoments_F_cumstd = np.zeros((nSubjects, len(reach_heights_simulated)))


"""
Running Simulations
"""

timestart = datetime.datetime.now()
print("Start time:-", timestart)

for sex_sim in range(0, len(sexs)):

    sex = sexs[sex_sim]

    # Center of Mass Location: % longitudinal position from joint center
    # Taken from: de Leva (1996). doi: 10.1016/0021-9290(95)00178-6 (Table 4)
    if sex == 'M':
        segment_com_loc = [
            39.50,  # Hand (wrt wrist joint centre; calculated as 79.00/2 because endpoint was to the palm)
            45.74,  # Forearm (wrt elbow joint centre; proximal)
            57.72,  # Upperarm (wrt shoulder joint centre; proximal)
            55.14,  # Torso (wrt mid-hip; distal)
            59.05,  # Thigh (wrt knee joint centre; distal)
            55.41   # Shank (wrt lateral malleolus; distal)
        ]
    elif sex == 'F':
        segment_com_loc = [
            37.37,  # Hand (wrt wrist joint centre; calculated as 74.74/2 because endpoint was to the palm)
            45.59,  # Forearm (wrt elbow joint centre; proximal)
            57.54,  # Upperarm (wrt shoulder joint centre; proximal)
            58.49,  # Torso (wrt mid-hip; distal)
            63.88,  # Thigh (wrt knee joint centre; distal)
            55.84   # Shank (wrt lateral malleolus; distal)
        ]

    # Mean Isometric Joint Strengths
    # Wrist F/E: Holzbaur et al. (2007). DOI: 10.1016/j.jbiomech.2006.11.013 (Table 1)
    # Wrist RD/UD (M): Delp et al. (1996). DOI: 10.1016/0021-9290(96)00029-2 (Table 1)
    # Wrist RD/UD (F): Multiplied values from males based on average sex difference from Plewa et al. (2016): http://dx.doi.org/10.1016/j.apergo.2015.09.005
    # Elbow / Shoulder: Kumar (2004). Muscle strength (Pages 241-243)
    # Trunk: Keller & Roy (2002). DOI: 10.1097/00024720-200208000-00009 (Figure 2)
    if sex == 'M':
        #strength_wrist_flx = 25.6
        #strength_wrist_ext = 14
        strength_wrist_rad = 11.0
        strength_wrist_uln = 9.1
        strength_elb_flx = 77
        strength_elb_ext = 46
        strength_shd_flx = 68
        strength_shd_ext = 93
        strength_tnk_flx = 123
        strength_tnk_ext = 218
    elif sex == 'F':
        #strength_wrist_flx = 10.7
        #strength_wrist_ext = 6.4
        strength_wrist_rad = 11.0 * (3/5)
        strength_wrist_uln = 9.1 * (3/5)
        strength_elb_flx = 41
        strength_elb_ext = 27
        strength_shd_flx = 45
        strength_shd_ext = 54
        strength_tnk_flx = 60
        strength_tnk_ext = 137

    for subject in range(0, nSubjects):

        print('--------')
        print("Subject #", subject + 1, sex)

        '''
        1. Determine Segment Properties
        - From proportionality constants
        - Given: Sex 'M' or 'F'
        '''

        # Return anthropometrics
        anthro = segments.segmentparameters(sex, average=evalaverage)  # Could sample nSubjects all at once instead of 1 at a time for speed
        segment_lengths = anthro[0]
        segment_mass = anthro[1]
        subject_mass = anthro[3]
        subject_height = anthro[4]

        if constrainRelative is True:
            armLength = anthro[0][0] + anthro[0][1] + anthro[0][2]
            constrainDistance = armLength * (constrainRelativePercent / 100)
        else:
            pass

        for reach_sim in range(0, len(reach_heights_simulated)):

            for joint_sim in range(0, len(joints_perturbed)):

                reach_height = reach_heights_simulated[reach_sim]
                joint = joints_perturbed[joint_sim]

                # print("Starting: ",
                #      " reach height ", reach_heights_simulated[reach_sim],
                #      " joint ", joints_perturbed[joint_sim])

                '''
                2. Defining Initial Positions (needed for optimization)
                '''
                theta_0 = [
                    0,      # Hand (pulling or right-side) segment angle
                    45,     # Forearm (pulling or right-side) segment angle
                    315,    # Upperarm (pulling or right-side) segment angle
                    315,    # Hand (non-pulling or left-side) segment angle
                    315,    # Forearm (non-pulling or left-side) segment angle
                    315,    # Upperarm (non-pulling or left-side) segment angle
                    80,     # Trunk segment angle
                    100,    # Thigh segment angle
                    80,     # Shank segment angle
                    ]

                if visualizeInitialguess is True:
                    positions = posture.jointpositions(theta_0, segment_lengths)
                    bodycom_x = posture.bodycenterofmass(positions[0], segment_mass, segment_com_loc)
                    visualize.visualizeposture(positions[0], positions[1], bodycom_x, time=0.2)
                else:
                    pass

                '''
                3. Defining Constraints
                - Joint angle limits
                - Horizontal Distance (change depending on individual stature or keep absolute like heights)
                - Vertical Distance (change depending on 4 heights)
                - Balance
                    - make sure COP is within the base of support
                '''

                # 3a. Joint Angles
                # Each constrained to their physiological range
                # Equality constraint means that the constraint function result is to be zero (= 0)
                # Inequality means that it is to be non-negative (>= 0)

                # Joint angle limits
                # Note: numbers here expressed/reported in different conventions for interpretation (see output data section)
                # Wrist: -25 <= theta[1]-theta[0] <= 40         -25 = radial deviation; 40 = ulnar deviation
                # Elbow: 210 <= theta[2]-theta[1] <= 360        210 = 30 degree elbow flexion; 360 = fully extended elbow
                # Shoulder: -360 <= theta[6]-theta[2] <= -150   -360 = arm flexed fwd 180; -150 = arm extended back 30 deg
                # Torso: -30 <= theta[7]-theta[6] <= 180          0 = straight trunk upright; 180 = trunk fully flexed
                # Knee: -180 <= theta[8]-theta[7] <= 0          -180 = fully flexed knee; 0 = fully extended knee
                # Non-pulling (left-side) arm constrained to pointing down: theta[3], theta[4], theta[5] = 270

                con1 = {'type': 'ineq', 'fun': lambda theta: theta[0]}  # straight hand
                con2 = {'type': 'ineq', 'fun': lambda theta: -theta[0]}  # straight hand
                con3 = {'type': 'ineq', 'fun': lambda theta: theta[1] - theta[0] + 25}  # radial deviation
                con4 = {'type': 'ineq', 'fun': lambda theta: -theta[1] + theta[0] + 40}  # ulnar deviation
                con5 = {'type': 'ineq', 'fun': lambda theta: theta[2] - theta[1] - 210}  # elbow min
                con6 = {'type': 'ineq', 'fun': lambda theta: -theta[2] + theta[1] + 360}  # elbow max
                con7 = {'type': 'ineq', 'fun': lambda theta: -theta[6] + theta[2] - 150}  # shoulder min
                con8 = {'type': 'ineq', 'fun': lambda theta: theta[6] - theta[2] + 360}  # shoulder max
                con9 = {'type': 'ineq', 'fun': lambda theta: -theta[7] + theta[6] + 180}  # torso min
                con10 = {'type': 'ineq', 'fun': lambda theta: theta[7] - theta[6] + 30}  # torso max
                con11 = {'type': 'ineq', 'fun': lambda theta: theta[8] - theta[7] + 180}  # knee min
                con12 = {'type': 'ineq', 'fun': lambda theta: -theta[8] + theta[7]}  # knee max
                con13 = {'type': 'eq', 'fun': lambda theta: theta[3] - 270}  # wrist (non-pulling or left-side)
                con14 = {'type': 'eq', 'fun': lambda theta: theta[4] - 270}  # elbow (non-pulling or left-side)
                con15 = {'type': 'eq', 'fun': lambda theta: theta[5] - 270}  # shoulder (non-pulling or left-side)

                # 3b. Ankle Location
                # Alter x position for horizontal distance and y position for vertical distance
                # Initial joint positions (m) relative to hand endpoint: x = 0.3929437746996767, y = -0.738880968349048
                # con16 and con17 refer to the horizontal distance of the ankle
                # constrainFoot = True: foot is at an exact horizontal distance away from hand (absolute or relative)
                # constrainFoot = False: foot is within 0-200 cm horizontal distance away from hand
                # con18 enforces that the reach height of the task is constrained

                if constrainFoot is True:
                    def horizontalconmin(theta):
                        bodyposition = posture.jointpositions(theta, segment_lengths)
                        return bodyposition[0][5] - constrainDistance

                    def horizontalconmax(theta):
                        bodyposition = posture.jointpositions(theta, segment_lengths)
                        return constrainDistance - bodyposition[0][5]

                    def verticalcon(theta):
                        bodyposition = posture.jointpositions(theta, segment_lengths)
                        return bodyposition[1][5] - anthro[2] + reach_height

                else:
                    def horizontalconmin(theta):
                        bodyposition = posture.jointpositions(theta, segment_lengths)
                        return bodyposition[0][5] - 0

                    def horizontalconmax(theta):
                        bodyposition = posture.jointpositions(theta, segment_lengths)
                        return 200 - bodyposition[0][5]

                    def verticalcon(theta):
                        bodyposition = posture.jointpositions(theta, segment_lengths)
                        return bodyposition[1][5] - anthro[2] + reach_height

                con16 = {'type': 'ineq', 'fun': horizontalconmin}
                con17 = {'type': 'ineq', 'fun': horizontalconmax}
                con18 = {'type': 'eq', 'fun': verticalcon}

                # 3c. Balance
                # COP must be inside base of support or else they fall over
                # Enforced as body's COM must be same horizontal distance from hand as the ankle

                def balanceconfwd(theta):
                    bodyposition = posture.jointpositions(theta, segment_lengths)
                    bodycom = posture.bodycenterofmass(bodyposition[0], segment_mass, segment_com_loc)
                    return bodycom - bodyposition[0][5]

                def balanceconbwk(theta):
                    bodyposition = posture.jointpositions(theta, segment_lengths)
                    bodycom = posture.bodycenterofmass(bodyposition[0], segment_mass, segment_com_loc)
                    return bodyposition[0][5] - bodycom

                con19 = {'type': 'ineq', 'fun': balanceconfwd}
                con20 = {'type': 'ineq', 'fun': balanceconbwk}

                # Defining array of constraints
                cons = [ con1, con2, con3, con4, con5, con6, con7, con8, con9, con10,
                con11, con12, con13, con14, con15, con16, con17, con18, con19, con20 ]

                '''
                4. Initial Optimization
                '''

                # Objective Function
                def objective(theta):

                    sum_squares = math.sqrt(sum(np.array(moments.inversedynamics1(theta, segment_mass, segment_lengths,
                                                                        segment_com_loc, load, loadangle))**2))
                    return sum_squares

                # Initial Minimal Solution (kinematics in degrees)
                minimize_sol = minimize(objective, theta_0, method='SLSQP', constraints=cons)
                sol_angles = minimize_sol.x
                # print("Initial solution angles (in degrees)", sol_angles)

                positions = posture.jointpositions(sol_angles, segment_lengths)
                bodycom_x = posture.bodycenterofmass(positions[0], segment_mass, segment_com_loc)
                #print(bodycom_x)

                if visualizeInitialsolution is True and joint_sim == 0:
                    if evalaverage is True:
                        if sex == 'M':
                            text = 'Male - ' + str(reach_height) + ' cm'
                        else:
                            text = 'Female - ' + str(reach_height) + ' cm'
                        saveplotname = 'Figures/Fig1-PosturePrediction/' + text + '.png'
                        visualize.visualizeposture(positions[0], positions[1], bodycom_x, time=0.5, title=text, plotfile=saveplotname)
                    else:
                        visualize.visualizeposture(positions[0], positions[1], bodycom_x, time=0.5)
                else:
                    pass

                '''
                5. Joint Perturbations
                - +/- 5 degree joint perturbations at specified joint
                - determines remaning joint angles subject to previous constraints
                - calculate moments from newly constrained joint position
                '''

                # Cycle though joint range of +/- perturbation degrees at specified joint
                for increase in range(-perturbation, perturbation+1):

                    # Joint Constraints
                    def wristperturbcon(joint_angles):
                        return (joint_angles[1] - (joint_angles[0] - 180)) \
                               - (sol_angles[1] - (sol_angles[0] - 180) + increase)

                    def elbowperturbcon(joint_angles):
                        return (joint_angles[1] * -1 + (joint_angles[2] - 180)) \
                               - ((sol_angles[1] * -1 + (sol_angles[2] - 180)) + increase)

                    def shoulderperturbcon(joint_angles):
                        old = joint_angles[2] - (joint_angles[6] + 180)
                        new = sol_angles[2] - (sol_angles[6] + 180)
                        diff = old - new
                        return diff - increase

                    # Setting constraints to lock horizontal location of ankle to same as initial optimized solution (optional)
                    def horizontalconmin_lock(theta):
                        bodyposition = posture.jointpositions(theta, segment_lengths)
                        return bodyposition[0][5] - bodycom_x

                    def horizontalconmax_lock(theta):
                        bodyposition = posture.jointpositions(theta, segment_lengths)
                        return bodycom_x - bodyposition[0][5]

                    con16_lock = {'type': 'ineq', 'fun': horizontalconmin_lock}
                    con17_lock = {'type': 'ineq', 'fun': horizontalconmax_lock}

                    if joint == 'Wrist':
                        angleCon = {'type': 'eq', 'fun': wristperturbcon}
                        consperturb = [ con1, con2, angleCon, con5, con6, con7, con8, con9, con10, con11, con12, con13, con14, con15, con16, con17, con18, con19, con20 ]

                    elif joint == 'Elbow':
                        angleCon = {'type': 'eq', 'fun': elbowperturbcon}
                        consperturb = [ con1, con2, con3, con4, angleCon, con7, con8, con9, con10, con11, con12, con13, con14, con15, con16, con17, con18, con19, con20 ]

                    elif joint == 'Shoulder':
                        angleCon = {'type': 'eq', 'fun': shoulderperturbcon}
                        consperturb = [ con1, con2, con3, con4, con5, con6, angleCon, con9, con10, con11, con12, con13, con14, con15, con16, con17, con18, con19, con20 ]

                    # Change from optimal solution
                    change_angle_solution = minimize(objective, sol_angles, method='SLSQP', constraints=consperturb)
                    change_angle_sol = change_angle_solution.x

                    # Reporting external moments (inverse dynamics here calculates internal moments)
                    idmoments = -1*np.array(moments.inversedynamics1(change_angle_sol, segment_mass, segment_lengths,
                                                                       segment_com_loc, load, loadangle))
                    idmoments[7:] = idmoments[7:] / 2  # Distributing knee and ankle moments across two sides
                    rw_moment = idmoments[0]
                    re_moment = idmoments[1]
                    rs_moment = idmoments[2]
                    lw_moment = idmoments[3]
                    le_moment = idmoments[4]
                    ls_moment = idmoments[5]
                    b_moment = idmoments[6]
                    k_moment = idmoments[7]
                    a_moment = idmoments[8]

                    positions = posture.jointpositions(change_angle_solution.x, segment_lengths)
                    bodycom_x = posture.bodycenterofmass(positions[0], segment_mass, segment_com_loc)
                    #print(bodycom_x)

                    if visualizechange is True:
                        title = 'Subject # ' + str(subject) + '; Reach Height: ' + str(reach_height) + 'cm; ' + joint + ' ' + str(increase)
                        visualize.visualizeposture(positions[0], positions[1], bodycom_x, time=0.2, title=title)
                    else:
                        pass

                    # Normalizing the joint moments to average joint strength
                    if rw_moment > 0:
                        #norm_rw_moment = (rw_moment / strength_wrist_flx) * 100
                        norm_rw_moment = (rw_moment / strength_wrist_rad) * 100
                    else:
                        #norm_rw_moment = (rw_moment / strength_wrist_ext) * 100
                        norm_rw_moment = (rw_moment / strength_wrist_uln) * 100

                    if re_moment > 0:
                        norm_re_moment = (re_moment / strength_elb_flx) * 100
                    else:
                        norm_re_moment = (re_moment / strength_elb_ext) * 100

                    if rs_moment > 0:
                        norm_rs_moment = (rs_moment / strength_shd_flx) * 100
                    else:
                        norm_rs_moment = (rs_moment / strength_shd_ext) * 100

                    if b_moment > 0:
                        norm_b_moment = (b_moment / strength_tnk_ext) * 100
                    else:
                        norm_b_moment = (b_moment / strength_tnk_flx) * 100


                    """
                    DATA OUTPUT
                    """
                    if sex == 'M':
                        subjectid = subject
                    elif sex == 'F':
                        subjectid = subject + nSubjects

                    # Error handling for failed optimizations
                    if change_angle_solution.success is False:
                        simulationdata = [
                            subjectid + 1,
                            sex,
                            subject_mass,
                            subject_height,
                            reach_height,
                            joint,     # Joint Perturbed
                            increase,  # Joint Perturbation Magnitude
                            None,  # Wrist Moment
                            None,  # Elbow Moment
                            None,  # Shoulder Moment
                            None,  # Back Moment
                            None,  # Normalized Wrist Moment
                            None,  # Normalized Elbow Moment
                            None,  # Normalized Shoulder Moment
                            None,  # Normalized Back moment
                            None,  # Wrist Angle
                            None,  # Elbow Angle
                            None,  # Shoulder Angle
                            None   # Trunk Angle
                        ]
                        failedsolutions += 1

                    # Successful optimizations (joint angles expressed in right-hand rule convention)
                    else:
                        simulationdata = [
                            subjectid + 1,
                            sex,
                            subject_mass,
                            subject_height,
                            reach_height,
                            joint,
                            increase,
                            rw_moment,
                            re_moment,
                            rs_moment,
                            b_moment,
                            norm_rw_moment,
                            norm_re_moment,
                            norm_rs_moment,
                            norm_b_moment,
                            (change_angle_sol[1] - change_angle_sol[0]),                # Wrist Angle (+ = ulnar deviation; - = radial deviation)
                            (change_angle_sol[2] - change_angle_sol[1] - 180),          # Elbow Angle (0 = fully flexed; 180 = fully extended)
                            (change_angle_sol[6] - change_angle_sol[2] + 180),          # Shoulder Angle (+ = extension; - = flexion)
                            (change_angle_sol[7] - change_angle_sol[6]),                # Trunk Angle (0 = upright trunk; + = flexion)
                        ]

                        simulationresults.append(simulationdata)


                    """
                    STORING DATA FOR CONVERGENCE CHECK
                    """

                    if (joint == 'Wrist') & (increase == 0):
                        if sex == 'M':
                            optimal_rwmoments_M[subject, reach_sim] = rw_moment
                            optimal_remoments_M[subject, reach_sim] = re_moment
                            optimal_rsmoments_M[subject, reach_sim] = rs_moment
                            optimal_bmoments_M[subject, reach_sim] = b_moment

                            optimal_rwmoments_M_cummean[subject, reach_sim] = np.mean(optimal_rwmoments_M[0:subject+1, reach_sim])
                            optimal_remoments_M_cummean[subject, reach_sim] = np.mean(optimal_remoments_M[0:subject+1, reach_sim])
                            optimal_rsmoments_M_cummean[subject, reach_sim] = np.mean(optimal_rsmoments_M[0:subject+1, reach_sim])
                            optimal_bmoments_M_cummean[subject, reach_sim] = np.mean(optimal_bmoments_M[0:subject+1, reach_sim])

                            if subject > 0:
                                optimal_rwmoments_M_cumstd[subject, reach_sim] = np.std(optimal_rwmoments_M[0:subject+1, reach_sim], ddof=1)
                                optimal_remoments_M_cumstd[subject, reach_sim] = np.std(optimal_remoments_M[0:subject+1, reach_sim], ddof=1)
                                optimal_rsmoments_M_cumstd[subject, reach_sim] = np.std(optimal_rsmoments_M[0:subject+1, reach_sim], ddof=1)
                                optimal_bmoments_M_cumstd[subject, reach_sim] = np.std(optimal_bmoments_M[0:subject+1, reach_sim], ddof=1)

                        elif sex == 'F':
                            optimal_rwmoments_F[subject, reach_sim] = rw_moment
                            optimal_remoments_F[subject, reach_sim] = re_moment
                            optimal_rsmoments_F[subject, reach_sim] = rs_moment
                            optimal_bmoments_F[subject, reach_sim] = b_moment

                            optimal_rwmoments_F_cummean[subject, reach_sim] = np.mean(optimal_rwmoments_F[0:subject+1, reach_sim])
                            optimal_remoments_F_cummean[subject, reach_sim] = np.mean(optimal_remoments_F[0:subject+1, reach_sim])
                            optimal_rsmoments_F_cummean[subject, reach_sim] = np.mean(optimal_rsmoments_F[0:subject+1, reach_sim])
                            optimal_bmoments_F_cummean[subject, reach_sim] = np.mean(optimal_bmoments_F[0:subject+1, reach_sim])

                            if subject > 0:
                                optimal_rwmoments_F_cumstd[subject, reach_sim] = np.std(optimal_rwmoments_F[0:subject+1, reach_sim], ddof=1)
                                optimal_remoments_F_cumstd[subject, reach_sim] = np.std(optimal_remoments_F[0:subject+1, reach_sim], ddof=1)
                                optimal_rsmoments_F_cumstd[subject, reach_sim] = np.std(optimal_rsmoments_F[0:subject+1, reach_sim], ddof=1)
                                optimal_bmoments_F_cumstd[subject, reach_sim] = np.std(optimal_bmoments_F[0:subject+1, reach_sim], ddof=1)


"""
SAVING DATA
"""
if evalaverage is True:
    imageset = ['Figures/Fig1-PosturePrediction/Male - 140 cm.png', 'Figures/Fig1-PosturePrediction/Female - 140 cm.png', 'Figures/Fig1-PosturePrediction/Male - 120 cm.png', 'Figures/Fig1-PosturePrediction/Female - 120 cm.png',
                'Figures/Fig1-PosturePrediction/Male - 100 cm.png', 'Figures/Fig1-PosturePrediction/Female - 100 cm.png', 'Figures/Fig1-PosturePrediction/Male - 80 cm.png', 'Figures/Fig1-PosturePrediction/Female - 80 cm.png']
    visualize.combine_images(columns=2, space=0, images=imageset)
    sys.exit()

simulationresults_pd = pd.DataFrame(simulationresults,
                                    columns=['Subject', 'Sex', 'Subject_Mass', 'Subject_Height',
                                             'Reach_Height', 'Joint', 'Joint_Pert',
                                             'RW_Moment', 'RE_Moment', 'RS_Moment', 'B_Moment',
                                             'Norm_RW_Moment', 'Norm_RE_Moment', 'Norm_RS_Moment', 'Norm_B_Moment',
                                             'Wrist_Angle', 'Elbow_Angle', 'Shoulder_Angle', 'Trunk_Angle'])
simulationresults_pd.to_csv('Data_Output/simulationresults.csv', encoding='utf-8')


"""
CHECK SAMPLING OF ANTHROPOMETRICS
"""
# Population Survey Data
popndata = pd.read_csv('samadult.csv', usecols=['AHEIGHT', 'SEX', 'AWEIGHTP'])
heightid = popndata[popndata['AHEIGHT'] > 85].index
popndata.drop(heightid, inplace=True)
weightid = popndata[popndata['AWEIGHTP'] > 800].index
popndata.drop(weightid, inplace=True)
popndata['AHEIGHT'] = popndata['AHEIGHT'] * 2.54
popndata['AWEIGHTP'] = popndata['AWEIGHTP'] * 0.453592
popndata_M = popndata[popndata['SEX'] == 1]
popndata_F = popndata[popndata['SEX'] == 2]

# Sampled Data
import matplotlib.pyplot as plt
male = simulationresults_pd[(simulationresults_pd.Sex == "M")]
female = simulationresults_pd[(simulationresults_pd.Sex == "F")]

# Plot
fig, axs = plt.subplots(2, sharex=True, sharey=True)
axs[0].scatter(x=popndata_M['AHEIGHT'], y=popndata_M['AWEIGHTP'], c='blue', alpha=0.05)
axs[0].scatter(x=popndata_F['AHEIGHT'], y=popndata_F['AWEIGHTP'], c='red', alpha=0.05)
axs[1].scatter(x=male['Subject_Height'], y=male['Subject_Mass'], c='white', edgecolor='blue')
axs[1].scatter(x=female['Subject_Height'], y=female['Subject_Mass'], c='white', edgecolor='red')
axs[0].set(title='Population Data', ylabel='Mass (kg)')
axs[1].set(title='Sampled Data', xlabel='Height (cm)', ylabel='Mass (kg)')
plt.savefig('Convergence/SamplingAnthropometrics')


"""
CONVERGENCE CHECK
"""
# Check change in running mean over last n simulations
nsimcheck = 50
compiledrunningmeans = np.array([optimal_rwmoments_M_cummean[-nsimcheck:, :], optimal_remoments_M_cummean[-nsimcheck:, :],
    optimal_rsmoments_M_cummean[-nsimcheck:, :], optimal_bmoments_M_cummean[-nsimcheck:, :],
    optimal_rwmoments_F_cummean[-nsimcheck:, :], optimal_remoments_F_cummean[-nsimcheck:, :],
    optimal_rsmoments_F_cummean[-nsimcheck:, :], optimal_bmoments_F_cummean[-nsimcheck:, :]])

maxrunningmean = np.amax(compiledrunningmeans, axis=1)
minrunningmean = np.amin(compiledrunningmeans, axis=1)
deltarunningmean = maxrunningmean-minrunningmean
np.set_printoptions(suppress=True)
print('Maximum change over the last ' + str(nsimcheck) + ' simulations in running mean of moment (Nm) across all joints and both sexes (rows) for all reaching heights (columns):')
print(deltarunningmean)

# Visual inspection of running mean and standard deviation
subject_idx = np.arange(1, nSubjects+1)

for sex_sim in range(0, len(sexs)):

    for reach_sim in range(0, len(reach_heights_simulated)):

        fig, axs = plt.subplots(4, 2, figsize=(15, 15))
        plot_title = sexs[sex_sim] + ': ' + str(reach_heights_simulated[reach_sim]) + ' cm'
        fig.suptitle(plot_title, fontweight='bold', fontsize=40.0)

        if sexs[sex_sim] == 'M':

            axs[0, 0].plot(subject_idx, optimal_rwmoments_M_cummean[:, reach_sim])
            axs[1, 0].plot(subject_idx, optimal_remoments_M_cummean[:, reach_sim])
            axs[2, 0].plot(subject_idx, optimal_rsmoments_M_cummean[:, reach_sim])
            axs[3, 0].plot(subject_idx, optimal_bmoments_M_cummean[:, reach_sim])
            axs[0, 0].set(title='Mean', ylabel='Wrist Moment (Nm)')
            axs[1, 0].set(ylabel='Elbow Moment (Nm)')
            axs[2, 0].set(ylabel='Shoulder Moment (Nm)')
            axs[3, 0].set(ylabel='Trunk Moment (Nm)', xlabel='# of Subjects')

            axs[0, 1].plot(subject_idx[1:], optimal_rwmoments_M_cumstd[1:, reach_sim])
            axs[1, 1].plot(subject_idx[1:], optimal_remoments_M_cumstd[1:, reach_sim])
            axs[2, 1].plot(subject_idx[1:], optimal_rsmoments_M_cumstd[1:, reach_sim])
            axs[3, 1].plot(subject_idx[1:], optimal_bmoments_M_cumstd[1:, reach_sim])
            axs[0, 1].set(title='Stdev')
            axs[3, 1].set(xlabel='# of Subjects')


        elif sex == 'F':

            axs[0, 0].plot(subject_idx, optimal_rwmoments_F_cummean[:, reach_sim])
            axs[1, 0].plot(subject_idx, optimal_remoments_F_cummean[:, reach_sim])
            axs[2, 0].plot(subject_idx, optimal_rsmoments_F_cummean[:, reach_sim])
            axs[3, 0].plot(subject_idx, optimal_bmoments_F_cummean[:, reach_sim])
            axs[0, 0].set(title='Mean', ylabel='Wrist Moment (Nm)')
            axs[1, 0].set(ylabel='Elbow Moment (Nm)')
            axs[2, 0].set(ylabel='Shoulder Moment (Nm)')
            axs[3, 0].set(ylabel='Trunk Moment (Nm)', xlabel='# of Subjects')

            axs[0, 1].set(title='Stdev')
            axs[0, 1].plot(subject_idx[1:], optimal_rwmoments_F_cumstd[1:, reach_sim])
            axs[1, 1].plot(subject_idx[1:], optimal_remoments_F_cumstd[1:, reach_sim])
            axs[2, 1].plot(subject_idx[1:], optimal_rsmoments_F_cumstd[1:, reach_sim])
            axs[3, 1].plot(subject_idx[1:], optimal_bmoments_F_cumstd[1:, reach_sim])
            axs[3, 1].set(xlabel='# of Subjects')


        plt.savefig('Convergence/' + 'convergence' + sexs[sex_sim] + str(reach_heights_simulated[reach_sim]))


print("Number of failed solutions = ", failedsolutions)
timeend = datetime.datetime.now()
print("Start time:-", timestart)
print("End time:-", timeend)
