# ==========================================================================================================================

# Spatio-Temporal Phase-Shifhting Profilometry algorithm

# This code is derived from the method detailed in the paper by Ri 2019 [1]; it uses the unwrapping algorithm of Herraez 2002 [2]
# and the phase-to-height relation derived by Maurel 2009 [3].

# References:
# [1] Ri S., Wang Q., Xia P., Tsuda H.: Spatiotemporal phase-shifting method for accurate phase analysis of fringe pattern.
# Journal of Optics 21(9), 095702 (2019)
# [2] Herraez M.A., Burton D.R., Lalor M.J., Gdeisat M.A.: Fast two-dimensional phase-unwrapping algorithm based on sorting 
# by reliability following a noncontinuous path. Applied optics 41(35), 7437–7444 (2002)
# [3] Maurel A., Cobelli P., Pagneux V., Petitjeans P.: Experimental and theoretical inspection of the phase-to-height 
# relation in fourier transform profilometry. Applied Optics 48(2), 380–392 (2009)

# ==========================================================================================================================

# In this code, a number N_shift = 3 of projected fringe patterns were used on both the reference plane and the object 
# of study (a rivulet flowing down a 50.4° incline at a flow rate Q=45mL/min), with a total number 2*N_shift = 6 images.

# All images should be stored in the path_im folder and the result of the algorithm is h.tif, the reconstructed height map 
# of the object.
# Here the "Ref" images correspond to the reference plane and the "Def" images to the object of study.
# The number following "Ref" or "Def" corresponds to the phase shift of the projected fringe pattern: 0 for a phase shift 
# of 0, 1 for a phase shift of 2*pi/3, 1 for a phase shift of 4*pi/3.

# The notations L, D and w0 are in agreement with the notations of Maurel [3] and correspond to the projecting distance (L),
# the distance between the camera and projector (D) and the fringe pattern frequency on the reference plane (w0).
# L and D were pre-calculated through a calibration using a solid wedge. The sampling period T_sampling corresponds to the
# period T as defined by Ri [1].

# ==========================================================================================================================

import numpy as np
import matplotlib.pyplot as plt
from unwrap import unwrap
from scipy import interpolate
import os
import cv2

#%% Parameter initialisation

N_shift = 3

T_sampling_float = 3143/100
T_sampling = round(T_sampling_float)

# Pixel ratio in mm/pix
R_pix = 100/3326.4

L, D = 804.74, 228.99
w0 = 2*np.pi/(R_pix*T_sampling_float)

# TO BE MODIFIED AFTER UPLOAD
path = '/Users/heliedemiramon/Desktop/TheseLadhyx/Papier Experiments in Fluids/Github/Rivulet 45 mL_min/'
path_im = path+'Exp/'

# Selecting and sorting the file names for all images

file_list = os.listdir(path_im)

for file in file_list:
    if not '.jpg' in file:
        file_list.remove(file)

file_list = sorted(file_list)

for t in range(2*N_shift):
    file_list[t] = path_im+file_list[t]

# Reading the images in greyscale using the cv2.imreading(file_name, 0) method; Img is defined as the list 
# containing all images as numpy.ndarrays

Img = [cv2.imread(file_list[k],0) for k in range(2*N_shift)]
(N_y,N_x) = Img[0].shape

# Normalization of the raw images; both A and B contain 2 arrays, one for the Def images, one for the Ref images

A = [1/N_shift*sum(Img[k+p*N_shift] for k in range(N_shift)) for p in range(2)]
B = [2/N_shift*(sum(Img[k+p*N_shift]*np.sin(2*np.pi*k/N_shift) for k in range(N_shift))**2+\
                    sum(Img[k+p*N_shift]*np.cos(2*np.pi*k/N_shift) for k in range(N_shift))**2)**0.5 for p in range(2)]

for n in range(2*N_shift):
    Img[n] = (Img[n]-A[n//N_shift])/B[n//N_shift]

# Down-sampling and interpolation calculations; all of these steps are detailed in Ri [1]. I_interpol contains 
# 2*N_shift*T_sampling arrays of size (N_y,N_x) corresponding to the down sampled and interpolated intensity profiles for
# both the object of study (p=0) and the reference plane (p=1).
# The resulting array Psi contains the overall phase for both the reference plane (psi_ref) and the object (psi)
# Both psi and psi_ref are then unwrapped using the 2d unwrapping algorithm of Herraez [2]

I_interpol = np.zeros((2, N_shift, T_sampling, N_y, N_x))
for p in range(2):
    for k in range(N_shift):
        for t in range(T_sampling):
            X = T_sampling*np.arange(N_x//T_sampling)+t
            Y = Img[N_shift*p+k][:,X]
            f_linear = interpolate.make_interp_spline(X, Y, k=1, axis=1)
            I_interpol[p,k,t] = f_linear(np.arange(N_x))

Psi_z = [sum(I_interpol[p,k,t]*np.exp(2*1j*np.pi*(t/T_sampling+k/N_shift)) \
             for k in range(N_shift) for t in range(T_sampling)) for p in range(2)]

Psi = -np.arctan2(np.real(Psi_z), np.imag(Psi_z))
[psi, psi_ref] = [unwrap(Psi[p]) for p in range(2)]

# Since the relevant information for 3d reconstructions is the phase difference between the object and the
# reference plane, there is no need to add a linear Moiré phase; one just to need to compoute the difference
# between psi and psi_ref; the height field h is then calculated with the phase-to-height relation of Maurel [3] 

delta_phi = psi - psi_ref

# The interpolation is erroneous on the first/last T_sampling columns (since there is only one point available for the
# linear polynomial interpolation); one can then simply crop the image on the left/right to avoid non-physical results on 
# the left/right edges of the phase difference map.
delta_phi = delta_phi[:,T_sampling:-T_sampling]

# Since the unwrapping algorithm defines both psi and psi_ref modulo a constant, one needs to define the zero on the
# phase difference map; in the following rivulet example, the height map of the object relative to the reference plane
# is 0 on the top left corner of the image, ie in a 10x10 pixels zone. The averaged value in this zone is subtracted from
# the overall phase difference map.
delta_phi -= np.mean(delta_phi[:10,:10])

# Calculation of the height map with the phase-to-height relation of Maurel [3] (neglecting the spatial parallax error,
# ie approximating (x',y')=(x,y) in equation 1.2 of the paper, since h/L<<1)
h = L*delta_phi/(delta_phi-w0*D)

plt.imshow(h, cmap='inferno')
plt.colorbar(shrink=0.7, label='h (mm)', location='bottom', pad=0.05)
plt.xticks([])
plt.yticks([])
plt.title('Reconstructed height map of the rivulet', pad=10)
plt.show()