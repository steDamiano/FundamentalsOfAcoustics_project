# FundamentalsOfAcoustics_project
Final project of Fundamentals of Acoustics Course @ Polimi (M. Sc. Music and Acoustic Engineering)

# Project aim
The goal of this project is to better analyze and understand the spectral behaviour of the Piano. It is well known that the vibration of a real string, with a finite mass and stiffness, departs from the harmonic behaviour, described by the equation:

![equation](https://latex.codecogs.com/gif.latex?f_%7Bn%7D%20%3D%20nf_%7B0%7D)

The stiff string actually produces harmonics with a frequency which is slightly sharp with respect to the ideal case: piano strings, specially in the low register, show a behaviour which is inharmonic to a different degree according to the piano register of the tone. The equation which describes the inharmonic progression of the partials is:

![equation](https://latex.codecogs.com/gif.latex?f_%7Bn%7D%20%3D%20nf_%7B0%7D%5Csqrt%7B1&plus;Bn%5E2%7D)

where B is the inharmonicity coefficient depending on string radius, linear density, tension, length and Young's modulus. 

The project started by taking piano recordings of two different instruments, a Grand Piano Steinway B2 and an Upright piano Yamaha U3, of which 4 notes of the low to middle register are sampled and analysed (C1, C2, C3, C4).
With this analysis the coefficient B is computed for each note and a comparison of the two instruments is made by showing the differences in the behaviour of the partials, using FFT, and of their decay, using STFT. 

# Algorithm used
The algorithm implemented for the analysis is proposed in:
```
Rauhala, Lehtonen, Välimäki
"Fast automatic inharmonicity estimation algorithm"
The Journal of the Acoustical Society of America, vol. 121, no. 5, pp. EL184-EL189, 2007.
```

The algorithm core loop is fully implemented, whereas the data is pre processed manually using Matlab FFT, without prabolic interpolation of the peaks. To make a comparison a naive algorithm is implemented too, using linear interpolation of the frequencies of the partials to get an estimate of the B coefficient.
Both algorithms are run over the samples of the Yamaha and the Steinway piano and results are shown with different graphs.

# Repository content
The repository has the Matlab files used to make data analysis together with samples used for the analysis. The Steinway samples are taken from Iowa University sitte at:
http://theremin.music.uiowa.edu/MISpiano.html<br>
Yamah U3 samples are recorded in a non-anechoic environment using two t.Bone SC140 Microphones placed 5cm over top of the piano, Left Mic over A2 and Right over A5. The recordings are stereo with a sample rate of 44100Hz.

The script "InharmonicityAnalysis.m" is the core of the project, implementing the Rauhala, Lehtonen Valimaki algorithm and plotting the results of the analysis of spectrum, partials and comparisons between ideal and real partials.
N.B: Audio Toolbox is required to run Secion II and hear the synthetic result of the sinusoidal spectrums build with ideal and real partials.
The script "PianoComparisons.m" concludes the analysis by comparing the inharmonicity coefficient of Upright and Grand Piano. It can be extended to include more Piano notes for a more accurate analysis.
