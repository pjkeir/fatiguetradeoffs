# fatiguetradeoffs

Team: Daanish Mulla, Nigel Majoni, Paul Tilley, Peter Keir

## Background
Upper body musculoskeletal disorders are the leading cause of lost time workplace claims in Canada.  These injuries are typically caused by repetitive, forceful tasks.  Lab-based studies investigating kinematic adaptations to fatiguing, repetitive work often find "subtle" joint angle changes (< 5 degrees).  We aimed to answer the following three questions:

1. What are the effects of anthropometric differences and task parameters on variability in joint demands?
2. Do meaningful changes in joint demands occur with "subtle" deviations in upper extremity joint angles?
3. How do changes in joint moments redistribute demands across the upper body?

## Approach
We used a probabilistic approach to generate a random sample of 1000 sagittal plane models based on population anthropometric variability.  For each model, we simulated a one-arm static pulling task at 50 N across four reach heights (80, 100, 120, and 140 cm) based on minimizing sum of joint moments.  To reflect "subtle" fatigue-induced changes, we altered upper extremity joints within +/- 5 degrees and re-solved for optimal postures.   

![me](https://github.com/pjkeir/fatiguetradeoffs/blob/main/Figures/gif/Posture%20Prediction.gif)
