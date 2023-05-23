# Light calibration with a leaf.

Download the executable files and the associated data in the [release page](https://github.com/Alienor134/OJIP-fit/releases/tag/preliminary_release)


1. Once the fluorescence curve recorded, save it as .csv file (separator: comma, decimal: point). If you don’t have a curve but want to try the app, use the OJIP_curve.csv provided in the .zip.  
2. Double click the executable file OJIP_fit.exe  
3. Open the address http://127.0.0.1:8050 on a web browser of your choice.  
4. Select the excitation wavelength used: it will display the associated sigma value.  
5. Drag-and-drop your .csv file (an animation is displayed while the file is loaded).  
6. Select the X and Y column names (you may need to click twice). An animation shows-up while the fit is being performed.  
7. The graph shows up with the 2 fitting methods. If you are not satisfied with the fit, play with the smoothing and logarithmic sub-sampling parameters until the beginning of the curve displays with a high point density.   
8. The tau value as well as the intensity values are displayed on the left. The error is expected to be a factor 2, which provides a reliable order of magnitude.   

![image](https://github.com/Alienor134/OJIP-fit/assets/20478886/388d8c00-0d01-4f13-b7a3-bd2be9eae28a)
