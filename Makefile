CC = mpic++
SRCS = *.cpp
INC = .
OPTS = -std=c++17 -Wall -Werror -O3

EXEC = nbody

all: clean compile

compile:
	$(CC) $(SRCS) $(OPTS) -I$(INC) -o $(EXEC)

clean:
	rm -f $(EXEC)