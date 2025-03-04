CC = mpicc
RM = /bin/rm
PROG = pingpong pingpong_i

all : ${PROG}

pingpong : pingpong.c
	${CC} -o $@ $^

pingpong_i : pingpong_i.c
	${CC} -o $@ $^

clean :
	${RM} -f ${PROG}