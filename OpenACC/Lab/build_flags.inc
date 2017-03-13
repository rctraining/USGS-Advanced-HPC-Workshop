FC = pgfortran
FCFLAGS_GPU_OMP = -acc -mp -fast -ta=tesla -Minfo=accel
FCFLAGS_GPU = -acc -fast -ta=tesla -Minfo=accel
FCFLAGS_SERIAL = -fast

