CC:=g++
exe:=main
obj:=main.o matrix.o

all:$(obj)
	$(CC)  -fopenmp -o $(exe) $(obj)
%.o:%.cpp 
	$(CC)  -fopenmp -c $^ -o $@ 

.PHONY:clean
clean:
	rm -rf $(obj) $(exe)

