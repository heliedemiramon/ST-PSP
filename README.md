# Spatio-Temporal Phase-Shifting Profilometry Algorithm

This code is derived from the method detailed in the paper by Ri 2019 [1]; it uses the unwrapping algorithm of Herraez 2002 [2] and the phase-to-height relation derived by Maurel 2009 [3].

## References

1. **Ri S., Wang Q., Xia P., Tsuda H.**: Spatiotemporal phase-shifting method for accurate phase analysis of fringe pattern. *Journal of Optics* **21**(9), 095702 (2019).

2. **Herraez M.A., Burton D.R., Lalor M.J., Gdeisat M.A.**: Fast two-dimensional phase-unwrapping algorithm based on sorting by reliability following a noncontinuous path. *Applied Optics* **41**(35), 7437–7444 (2002).

3. **Maurel A., Cobelli P., Pagneux V., Petitjeans P.**: Experimental and theoretical inspection of the phase-to-height relation in Fourier transform profilometry. *Applied Optics* **48**(2), 380–392 (2009).

---

## Implementation Details

- **Number of Fringe Patterns:**  
  A total of **N_shift = 3** projected fringe patterns were used on both the reference plane and the object of study. This results in a total of **2 × N_shift = 6 images**.

- **Image Storage and Output:**  
  - All images should be stored in the `path_im` folder.  
  - The result of the algorithm is a file named `h.tif`, which contains the reconstructed height map of the object.

- **Image Naming Convention:**  
  - **"Ref" images**: Correspond to the reference plane.  
  - **"Def" images**: Correspond to the object of study.  
  - The number following "Ref" or "Def" represents the phase shift of the projected fringe pattern:  
    - `0` for a phase shift of 0.  
    - `1` for a phase shift of 2*pi/3.  
    - `2` for a phase shift of 4*pi/3.

---

## Notations

The notations L, D, and w0 are consistent with those in Maurel [3]:  
- **L**: Projecting distance. 
- **D**: Distance between the camera and the projector.  
- **w0**: Fringe pattern frequency on the reference plane.

The sampling period \(**T_sampling**\) corresponds to the period **T** as defined by Ri [1].

---

## Calibration Details

Values for **L** and **D** were pre-calculated through a calibration process using a solid wedge. Detailed information about the calibration wedge is available in the `Calibration wedge` folder.
