# CC=/usr/local/gcc-5.4/bin/gcc
# CXX=/usr/local/gcc-5.4/bin/g++

CC=/usr/bin/gcc
CXX=/usr/bin/g++

SUB_DIR=$(shell ls -d */)

TARGET=$(SUB_DIR)

all:
	for dir in $(TARGET); \
		do \
		$(MAKE) -C $$dir CC=$(CC) CXX=$(CXX); \
		done

.PHONY: clean save test

clean:
	for dir in $(TARGET); \
		do \
			$(MAKE) clean -C $$dir; \
		done

# save the result in *.txt.bak to judge the correctness
save:
	for dir in $(SUB_DIR); \
		do \
			$(MAKE) save -C $$dir; \
		done

# test performance
test:
	for dir in $(SUB_DIR); \
		do \
			$(MAKE) test -C $$dir; \
		done
