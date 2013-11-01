function make
c = computer
switch c
    case 'i386-apple-darwin8.11.1'
        mex  -largeArrayDims -I. ...
                 fennelmex.c -D__thread= -DLINUX
    case 'MACI64'
        mex -O -largeArrayDims -I. ...
                 fennelmex.c -D__thread= -DLINUX
    case 'GLNXA64'
        mex -O -largeArrayDims -I. ...
                fennelmex.c -DLINUX CFLAGS="\$CFLAGS -std=c99" ...
                -D_POSIX_C_SOURCE=200809
    case 'GLNX32'
        mex -O -largeArrayDims -I. ...
                fennelmex.c -DLINUX
end
