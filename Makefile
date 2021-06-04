CC:= mpic++
exe:=main
obj:=main.o 

all:$(obj)
	$(CC) -o $(exe) $(obj)
%.o:%.cpp 
	$(CC) -c $^ -o $@ 

.PHONY:clean
clean:
	rm -rf $(obj) $(exe)

