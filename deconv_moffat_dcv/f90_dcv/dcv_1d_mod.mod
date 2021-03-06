V8 0 0 3 0
MODULE DCV_1D_MOD,2 0
FILE 0,dcv_1D_mod.f90
USE DCV_1D_DATA 2
REF IDIM(IDIM@DCV_1D_DATA),2
PROC DCV_GMARK,1,8,0: 8,0,0,0,0,262144,1,4194304
VAR ITER,3,0,,: 1,0,3,0,0,267,0,0
ENDPROC
PROC DCV_GGAUSS,1,8,0: 8,0,0,0,0,262144,1,4194304
VAR ITER,3,0,,: 1,0,3,0,0,267,0,0
ENDPROC
PROC DCV_TIKHONOV,1,8,0: 8,0,0,0,0,262144,1,4194304
VAR ITER,3,0,,: 1,0,3,0,0,267,0,0
ENDPROC
PROC DBRENT,7,8,0: 2,0,4,0,0,393345,2,4194304
VAR AX,3,0,,: 2,0,4,0,0,11,0,0
VAR BX,3,0,,: 2,0,4,0,0,259,0,0
VAR CX,3,0,,: 2,0,4,0,0,11,0,0
PROC FF,1,3,0: 2,0,4,0,0,131075,0,0
AAINFO 0 0 0 2 1
ENDPROC
PROC DFF,1,3,0: 2,0,4,0,0,131075,0,0
AAINFO 0 0 0 2 1
ENDPROC
VAR TOL,3,0,,: 2,0,4,0,0,259,0,0
VAR XMIN,3,0,,: 2,0,4,0,0,131,0,0
ENDPROC
PROC DCV_SQRTRF,2,8,0: 8,0,0,0,0,262144,1,4194304
VAR SS,3,0,,: 2,0,4,0,0,259,0,0
VAR ITER,3,0,,: 1,0,3,0,0,267,0,0
ENDPROC
PROC DCV_CONV_HH,4,8,0: 8,0,0,0,0,262144,2,0
VAR XX,3,0,,: 2,0,4,0,2397,8388867,0 (1,1,512: 1,512,1,1),0
VAR HH_RE,3,0,,: 2,0,4,0,2397,8388867,0 (1,1,512: 1,512,1,1),0
VAR HH_IM,3,0,,: 2,0,4,0,2397,8388867,0 (1,1,512: 1,512,1,1),0
VAR WW,3,0,,: 2,0,4,0,2397,8388747,0 (1,1,512: 1,512,1,1),0
ENDPROC
PROC MNBRAK,7,8,0: 8,0,0,0,0,393216,2,4194304
VAR AX,3,0,,: 2,0,4,0,0,395,0,0
VAR BX,3,0,,: 2,0,4,0,0,395,0,0
VAR CX,3,0,,: 2,0,4,0,0,395,0,0
VAR FA,3,0,,: 2,0,4,0,0,387,0,0
VAR FB,3,0,,: 2,0,4,0,0,387,0,0
VAR FC,3,0,,: 2,0,4,0,0,387,0,0
PROC FUNC,1,3,0: 2,0,4,0,0,131075,0,0
AAINFO 0 0 0 2 1
ENDPROC
ENDPROC
PROC DCV_MEM,1,8,0: 8,0,0,0,0,262144,1,4194304
VAR ITER,3,0,,: 1,0,3,0,0,267,0,0
ENDPROC
PROC F1DIM,1,8,0: 2,0,4,0,0,262281,2,4194304
VAR X,3,0,,: 2,0,4,0,0,259,0,0
ENDPROC
PROC DF1DIM,1,8,0: 2,0,4,0,0,262537,2,4194304
VAR X,3,0,,: 2,0,4,0,0,259,0,0
ENDPROC
PROC DCV_SPDIV,0,8,0: 8,0,0,0,0,262144,1,4194304
ENDPROC
PROC LINMIN,4,8,0: 8,0,0,0,0,393216,2,4194304
VAR P,3,0,,: 2,0,4,0,2956,8388995,0 (1,2,0: 1,2955,1,2),0
VAR XI,3,0,,: 2,0,4,0,2956,8388995,0 (1,2,0: 1,2955,1,2),0
VAR N,3,0,,: 1,0,3,0,0,259,0,0
VAR FRET,3,0,,: 2,0,4,0,0,131,0,0
ENDPROC
PROC DCV_WIENER,0,8,0: 8,0,0,0,0,262144,1,4194304
ENDPROC
PROC NORM_L2,2,8,0: 2,0,4,0,0,393601,2,4194304
VAR P,3,0,,: 2,0,4,0,2911,8388867,0 (1,2,0: 1,2910,1,2),0
VAR N,3,0,,: 1,0,3,0,0,259,0,0
ENDPROC
PROC FRPRMN,5,8,0: 8,0,0,0,0,393216,1,4194304
VAR P,3,0,,: 2,0,4,0,2604,267,0 (1,2,0: 1,2603,1,2),0
VAR N,3,0,,: 1,0,3,0,0,267,0,0
VAR FTOL,3,0,,: 2,0,4,0,0,259,0,0
VAR ITER,3,0,,: 1,0,3,0,0,459,0,0
VAR FRET,3,0,,: 2,0,4,0,0,267,0,0
ENDPROC
PROC DCV_DISPLAY_INPUT,0,8,0: 8,0,0,0,0,262144,1,4194304
ENDPROC
PROC DCV_PLOT2,4,8,0: 8,0,0,0,0,262144,1,4194304
VAR Y1,3,0,,: 2,0,4,0,241,8388619,0 (1,2,0: 1,240,1,2),0
VAR Y2,3,0,,: 2,0,4,0,241,8388619,0 (1,2,0: 1,240,1,2),0
VAR NN,3,0,,: 1,0,3,0,0,267,0,0
VAR TITLE,3,0,,: 3,0,10,40,0,16777227,0,0
ENDPROC
PROC DCV_PLOT1_LOG,3,8,0: 8,0,0,0,0,393216,1,4194304
VAR YY,3,0,,: 2,0,4,0,160,8388619,0 (1,2,0: 1,159,1,2),0
VAR NN,3,0,,: 1,0,3,0,0,259,0,0
VAR TITLE,3,0,,: 3,0,10,40,0,16777227,0,0
ENDPROC
PROC DCV_CHECK_GRAD,4,8,0: 8,0,0,0,0,393216,2,4194304
VAR X1,3,0,,: 2,0,4,0,2473,8388875,0 (1,1,512: 1,512,1,1),0
VAR X2,3,0,,: 2,0,4,0,2473,8388739,0 (1,1,512: 1,512,1,1),0
VAR DX,3,0,,: 2,0,4,0,2473,8388875,0 (1,1,512: 1,512,1,1),0
VAR NN,3,0,,: 1,0,3,0,0,259,0,0
ENDPROC
PROC DCV_PLOT1,3,8,0: 8,0,0,0,0,393216,1,4194304
VAR YY,3,0,,: 2,0,4,0,96,8388619,0 (1,2,0: 1,95,1,2),0
VAR NN,3,0,,: 1,0,3,0,0,267,0,0
VAR TITLE,3,0,,: 3,0,10,40,0,16777227,0,0
ENDPROC
PROC DFUNC_BANANA,2,8,0: 8,0,0,0,0,262144,2,4194304
VAR X,3,0,,: 2,0,4,0,1800,8388867,0 (1,1,2: 1,2,1,1),0
VAR DX,3,0,,: 2,0,4,0,1800,8388739,0 (1,1,2: 1,2,1,1),0
ENDPROC
PROC FUNC_BANANA,1,8,0: 2,0,4,0,0,262273,2,4194304
VAR X,3,0,,: 2,0,4,0,1764,8388867,0 (1,1,2: 1,2,1,1),0
ENDPROC
PROC DFUNC_MEM,2,8,0: 8,0,0,0,0,393216,2,4194304
VAR XX,3,0,,: 2,0,4,0,846,8389003,0 (1,1,512: 1,512,1,1),0
VAR DX,3,0,,: 2,0,4,0,846,8388747,0 (1,1,512: 1,512,1,1),0
ENDPROC
PROC FUNC_MEM,1,8,0: 2,0,4,0,0,393345,2,4194304
VAR XX,3,0,,: 2,0,4,0,774,8389003,0 (1,1,512: 1,512,1,1),0
ENDPROC
PROC DFUNC_TIKHO,2,8,0: 8,0,0,0,0,393216,2,4194304
VAR XX,3,0,,: 2,0,4,0,1094,8389003,0 (1,1,512: 1,512,1,1),0
VAR DX,3,0,,: 2,0,4,0,1094,8388747,0 (1,1,512: 1,512,1,1),0
ENDPROC
PROC FUNC_TIKHO,1,8,0: 2,0,4,0,0,393345,2,4194304
VAR XX,3,0,,: 2,0,4,0,1036,8389003,0 (1,1,512: 1,512,1,1),0
ENDPROC
PROC DFUNC_GMARK,2,8,0: 8,0,0,0,0,393216,2,4194304
VAR XX,3,0,,: 2,0,4,0,517,8389003,0 (1,1,512: 1,512,1,1),0
VAR DX,3,0,,: 2,0,4,0,517,8388747,0 (1,1,512: 1,512,1,1),0
ENDPROC
PROC FUNC_GMARK,1,8,0: 2,0,4,0,0,393345,2,4194304
VAR XX,3,0,,: 2,0,4,0,433,8389003,0 (1,1,512: 1,512,1,1),0
ENDPROC
PROC DFUNC_GGAUSS,2,8,0: 8,0,0,0,0,393216,2,4194304
VAR XX,3,0,,: 2,0,4,0,1329,8389003,0 (1,1,512: 1,512,1,1),0
VAR DX,3,0,,: 2,0,4,0,1329,8388747,0 (1,1,512: 1,512,1,1),0
ENDPROC
PROC DFUNC_SQRTRF,2,8,0: 8,0,0,0,0,393216,2,0
VAR XX,3,0,,: 2,0,4,0,1593,8389003,0 (1,1,512: 1,512,1,1),0
VAR DX,3,0,,: 2,0,4,0,1593,8388747,0 (1,1,512: 1,512,1,1),0
ENDPROC
PROC FUNC_GGAUSS,1,8,0: 2,0,4,0,0,393345,2,4194304
VAR XX,3,0,,: 2,0,4,0,1269,8388747,0 (1,1,512: 1,512,1,1),0
ENDPROC
PROC DFUNC,2,8,0: 8,0,0,0,0,393216,1,4194304
VAR X,3,0,,: 2,0,4,0,1917,8388619,0 (1,1,512: 1,512,1,1),0
VAR DX,3,0,,: 2,0,4,0,1917,8388619,0 (1,1,512: 1,512,1,1),0
ENDPROC
PROC FUNC,1,8,0: 2,0,4,0,0,393345,1,4194304
VAR X,3,0,,: 2,0,4,0,1853,8388619,0 (1,1,512: 1,512,1,1),0
ENDPROC
PROC FUNC_SQRTRF,1,8,0: 2,0,4,0,0,393345,2,0
VAR XX,3,0,,: 2,0,4,0,1531,8389003,0 (1,1,512: 1,512,1,1),0
ENDPROC
PROC DCV_BANANA,0,8,0: 8,0,0,0,0,262144,1,4194304
ENDPROC
END
