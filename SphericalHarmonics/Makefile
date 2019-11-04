CC	= g++
LD	= g++

#CCFLAGS	= -g -O0 `gsl-config --cflags` `root-config --cflags` -D_FINISH_DEBUG_
CCFLAGS	= -g -O0 `gsl-config --cflags` `root-config --cflags` 
LDFLAGS	= -g -O0 -lgsl `gsl-config --libs` `root-config --libs` 

testylmgsl: testylmgsl.o ylm.o mmatrix.o
	$(LD) -o $@ $^ $(LDFLAGS)

test: test.o ylm.o mmatrix.o
	$(LD) -o $@ $^ $(LDFLAGS)

%.o: %.cc
	$(CC) $^ -o $@ $(CCFLAGS) -c

clean:
	rm -f testylm *.o

testgslm: testgslm.o	
	$(LD) -o $@ $^ $(LDFLAGS)

testcg: testcg.o mmatrix.o
	$(LD) -o $@ $^ $(LDFLAGS)

testylm: testylm.o ylm.o
	$(LD) -o $@ $^ $(LDFLAGS)

CorrFctnDirectYlm.o: CorrFctnDirectYlm.cxx CorrFctnDirectYlm.h
	$(CC) $< -o $@ $(CCFLAGS) -c

testclass.o: testclass.cxx
	$(CC) $^ -o $@ $(CCFLAGS) -c

testclass: testclass.o CorrFctnDirectYlm.o ylm.o
	$(LD) -o $@ $^ $(LDFLAGS)

testclassread.o: testclassread.cxx
	$(CC) $^ -o $@ $(CCFLAGS) -c

testclassread: testclassread.o CorrFctnDirectYlm.o ylm.o
	$(LD) -o $@ $^ $(LDFLAGS)

testinteg.o: testinteg.cxx
	$(CC) $^ -o $@ $(CCFLAGS) -c

testinteg: testinteg.o CorrFctnDirectYlm.o ylm.o
	$(LD) -o $@ $^ $(LDFLAGS)

testcfcalc.o: testcfcalc.cxx
	$(CC) $^ -o $@ $(CCFLAGS) -c

testcfcalc: testcfcalc.o ylm.o
	$(LD) -o $@ $^ $(LDFLAGS)

integpion.o: integpion.cxx
	$(CC) $^ -o $@ $(CCFLAGS) -c

integpion: integpion.o ylm.o CorrFctnDirectYlm.o
	$(LD) -o $@ $^ $(LDFLAGS)

