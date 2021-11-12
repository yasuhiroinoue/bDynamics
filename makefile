# Makefile

TARGET = ./V_BDver03.another_diff.out
#TARGET = ./V_BDver03.sampling_rep.out
# TARGET = ./V_BDver03.sampling_rep_old_cluster.out
#TARGET = ./V_BDver03.nonpoly_rep.out

# --------For vt-opteron cluster---------
# CC = pgCC
# C_OPT = -fastsse -Msmart -Mvect=prefetch -Mipa=fast,inline -O3 -tp amd64
# CC = mpiCC
# CC = mpicxx
# C_OPT = -fastsse -Msmart -Mvect=prefetch -Mipa=fast,inline -O3 -tp amd64 -D=USE_MPI
# --------For old cluster option---------
# C_OPT = -fastsse -Msmart -Mvect=prefetch -Mipa=fast,inline -O3 -tp k8-64e

# --------For i7 cluster---------
# CC=icpc
# C_OPT=-fast -O3
#CC = mpicxx
#C_OPT=-ipo -O3 -D=USE_MPI

# --------For general cluster------------
CC = mpic++
C_OPT=-O3 -D=USE_MPI
#--------------------------------------

# --------For debug option---------
#C_OPT = 

# --------Object files-------------
OBJS = main.o _function.o mymath.o _arp2_3.o _class.o _pv_out.o _variable.o _output.o _output_class.o _mympi.o _wall_boundary.o _periodic.o

.SUFFIXES: .cpp .o

.cpp.o:
	$(CC) -c $(C_OPT) $<

.PHONY: all
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(C_OPT) -o $(TARGET) $(OBJS)
	@echo update change log!!
#	rm -f *.o
#	rm -f *.oo


$(OBJS): makefile *.h

.PHONY: clean
clean:
	@rm -f *.o 
	@rm -f *.oo
	@rm -f *.BAK

backupdirname = backup/${shell date +%m%d_%H%M}
backupdirname2 = ~/backup_code/${shell date +%Y%m}
.PHONY: backup
backup:
	@ mkdir -p $(backupdirname)/ &&mkdir -p $(backupdirname2)/ && cp *.h *.cpp makefile $(backupdirname)/ \
	&& tar czf $(backupdirname).tar.gz $(backupdirname)/ \
	&& cp  $(backupdirname).tar.gz $(backupdirname2)/\
	&&rm -fr $(backupdirname)/

