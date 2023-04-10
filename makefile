#main: main.c initSetting.c generateCACode.c acquisition.c tool.c track.c findPreamble.c satpos.c common.c cart2geo.c postProceeding.c
#	gcc main.c initSetting.c generateCACode.c tool.c acquisition.c track.c findPreamble.c satpos.c common.c cart2geo.c postProceeding.c -lgsl -lgslcblas -lm -o main -g -Wall


CC := gcc
CFLAGS := -lgsl -lgslcblas -lm -g -Wall
TARGET=GNSS_Receiver
OBJS=$(patsubst %.c, %.o, $(wildcard *.c))
.PHONY: clean

all: ${TARGET}

%.i: %.c
	$(info Preprocess: build *.i)
	${CC} ${CFLAGS} -E $< -o $@

%.s: %.c
	$(info Compile: build *.s)
	${CC} ${CFLAGS}-S $< -o $@

%.o: %.s
	$(info Asseble: build *.o)
	${CC} ${CFLAGS}-c $< -o $@

${TARGET}: ${OBJS} 
	$(info Link: build GNSS_Receiver)
	${CC} $^  ${CFLAGS} -o $@

clean:
	rm *.o


#afterTrackedMain: afterTrackedMain.c initSetting.c generateCACode.c acquisition.c tool.c track.c findPreamble.c satpos.c common.c cart2geo.c
#	gcc afterTrackedMain.c initSetting.c generateCACode.c tool.c acquisition.c track.c findPreamble.c satpos.c common.c cart2geo.c -lgsl -lgslcblas -lm -o main -g -Wall
