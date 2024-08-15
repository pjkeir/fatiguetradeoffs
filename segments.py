# Function to randomly sample anthropometrics to determine segment parameters for model creation
# Input: sex and population anthropometric data (see lines 11-13)
# Output: mass, height, and segment parameters (segment mass, segment lengths, COM location, inertial characteristics)
# Assumptions: proportionality constants based on height and mass determine segment parameters

import numpy as np
import pandas as pd

def segmentparameters(sex, average=False):

    # Population anthropometrics taken from the National Health Interview Survey (NHIS) on CDC website
    # Data available and downloaded from here: https://www.cdc.gov/nchs/nhis/nhis_2017_data_release.htm
    df = pd.read_csv('samadult.csv', usecols=['AHEIGHT', 'SEX', 'AWEIGHTP'])

    # Subset data based on sex (M or F)
    if sex == 'M':
        indexnames = df[df['SEX'] == 2].index
        df.drop(indexnames, inplace=True)
    else:
        indexnames = df[df['SEX'] == 1].index
        df.drop(indexnames, inplace=True)

    # Removing outliers / errors in datal
    heightid = df[df['AHEIGHT'] > 85].index
    df.drop(heightid, inplace=True)
    weightid = df[df['AWEIGHTP'] > 800].index
    df.drop(weightid, inplace=True)

    # Subset height and weight arrays of data
    height = df[['AHEIGHT']].to_numpy()
    weight = df[['AWEIGHTP']].to_numpy()
    height_weight = df[['AHEIGHT', 'AWEIGHTP']].to_numpy()

    # Determine means and covariance matrix of height and weight and creating random individual
    meanheight = np.mean(height)
    meanweight = np.mean(weight)
    cov = np.cov(height_weight.T)
    sub = np.random.multivariate_normal([meanheight, meanweight], cov)
    if average is True:
        sub[0] = meanheight
        sub[1] = meanweight

    # Converting to metric units (inches --> cm; lbs --> kg)
    subject_height = sub[0] * 2.54
    subject_mass = sub[1] * 0.453592
    print('Subject height: ', subject_height)
    print('Subject mass: ', subject_mass)

    # Proportionality constants for segment lengths based on body height
    # Taken from: Fromuth & Parkinson (2008). doi: 10.1115/DETC2008-50091 (Figure 5; 50th percentile)
    propotion_const_male = [
        .108 / 2,  # Fingertips to wrist (make sure location of load going through center of hand)
        .149,  # Wrist to elbow
        .191,  # Elbow to shoulder
        .290,  # Shoulder to low back
        .245,  # Low back to knee
        .246   # Knee to ankle
    ]

    propotion_const_female = [
        .110 / 2,  # Fingertips to wrist (make sure location of load going through center of hand)
        .153,  # Wrist to elbow
        .194,  # Elbow to shoulder
        .294,  # Shoulder to low back
        .242,  # Low back to knee
        .247   # Knee to ankle
    ]

    # Ankle proportions for reach heights
    if sex == 'M':
        ankle_height = 0.039 * subject_height
    elif sex == 'F':
        ankle_height = 0.038 * subject_height

    # Determine segment lengths from subject height
    segment_lengths = []

    if sex == 'M':
        for segment in range(0, len(propotion_const_male)):
            segment_lengths.append(subject_height * propotion_const_male[segment])
    elif sex == 'F':
        for segment in range(0, len(propotion_const_female)):
            segment_lengths.append(subject_height * propotion_const_female[segment])

    # Proportionality constants for segment masses based on body mass
    # Taken from: de Leva (1996). doi: 10.1016/0021-9290(95)00178-6
    mass_relation_male = [
        0.0061,  # Hand
        0.0162,  # Forearm
        0.0271,  # Upper arm
        0.4346,  # Torso
        0.1416 * 2,  # Thigh
        0.0433 * 2,  # Shank
    ]

    mass_relation_female = [
        0.0056,  # Hand
        0.0138,  # Forearm
        0.0255,  # Upper arm
        0.4257,  # Torso
        0.1478 * 2,  # Thigh
        0.0481 * 2,  # Shank
    ]

    # Determine segment masses from subject mass
    segment_mass = []

    if sex == 'M':
        for segment_ in range(0, len(mass_relation_male)):
            segment_mass.append(subject_mass * (mass_relation_male[segment_]))
    elif sex == 'F':
        for segment_ in range(0, len(mass_relation_female)):
            segment_mass.append(subject_mass * (mass_relation_female[segment_]))

    return[segment_lengths, segment_mass, ankle_height, subject_mass, subject_height]
