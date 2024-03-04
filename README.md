# Musical-Dissonance-Neural-Coding

Representation of consonance/dissonance of music intervals and chords by the midbrain neurons (i.e., inferior colliculus, IC). A modeling work. 

The initial results was shown in [Carney et al. ARO 2019 poster](https://github.com/TongShan4869/Musical-Dissonance-Neural-Coding/blob/main/ARO_2019_Carney_pitch_Final.pdf)



## Stimuli

**Experiment 1**

Dyads: root note A3 (f0=220 Hz) combining with the upper f0 in half semitone increment to create different intervals within the range of an octave - unison to octave.

**Experiment 2**

Dyads: root note A3 (f0=220 Hz) combining with the upper f0 in semitone increment to create different intervals within the range of an octave - unison to octave.

Triads: same root note A3 combining with uppter notes that form the primary chords in a western music within a key (power chord, major triad, minor triad, suspeded 4th, suspended 2nd, Italian 6th, diminished triad, augmented triad, Viennese trichord).

## Behaviroal data

**Experiment 1**

Behavioral consonance rating for 24 music intervals from Tufts et al. (2005). n=4.

**Experiment 2**

Behavioral consonance rating for 12 music intervals and 9 chords from Bowling et al. (2018). n=30.

## Modeled neurophysiology data

**SFIE model of Inferior Colliculus (Nelson & Carney, 2007)** was used to derive the IC neuron activity. Metrics to characterize consonance includes: 

1) standard deviation of mean neuron firing rate across characteristic frequencies (CF)
2) the standard deviation of neuron firing rate variation across CFs

The theory is based on the neural fluctuation mechanism of IC neurons (Carney, Li & McDonough, 2015). 

We then correlate the metrics with the behavioral data. 

![Consonance_correlation](https://github.com/TongShan4869/Musical-Dissonance-Neural-Coding/assets/51421789/9eef3770-c02c-462c-90e1-bb07d15127c0)

(Pearson's r=0.78)

![consonance_Bowling_chords](https://github.com/TongShan4869/Musical-Dissonance-Neural-Coding/assets/51421789/e01890a4-6eae-46f4-8520-5c3234991733)

(Pearson's r=0.94 for BE and r=0.85 for BS cell)

## Reference

Tufts, J. B. , Molis, M. R. , and Leek, M. R. (2005). “Perception of dissonance by people with normal hearing and sensorineural hearing loss,” J. Acoust. Soc. Am. 118, 955–967. https://doi.org/10.1121/1.1942347

Bowling, Daniel L., Dale Purves, and Kamraan Z. Gill. "Vocal similarity predicts the relative attraction of musical chords." Proceedings of the National Academy of Sciences 115.1 (2018): 216-221. 
https://doi.org/10.1073/pnas.171320611

Nelson PC, Carney LH (2007) Neural rate and timing cues for detection and discrimination of amplitude-modulated tones in the awake rabbit inferior colliculus. J Neurophysiol 97:522-539

Carney, Laurel H., Tianhao Li, and Joyce M. McDonough. "Speech coding in the brain: representation of vowel formants by midbrain neurons tuned to sound fluctuations." Eneuro 2.4 (2015).

