# Fast-High-Dimensional-Bilateral-and-Nonlocal-Means-Filtering
P. Nair, and K. N. Chaudhury, "Fast High-Dimensional Bilateral and Non-local Means Filtering," , IEEE Transactions on Image Processing, 28(3), pp.1470-1481.

IMPORTANT POINTS

Main files in the code:
1) bilateral_approx.m - For bilateral filtering of high-dimensional images.
2) nlm_approx.m - For non-local means filtering of high-dimensional images. 
3) bilateralHSdenoise.m - For hyperspectal filtering. 

If you need to use O(1) Gaussian filtering, you have to create mex file for young.cpp using the command
mex young.cpp 
