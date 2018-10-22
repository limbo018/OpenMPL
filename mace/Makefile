OBJ=mace

MACE: $(OBJ)

$(OBJ): mace.c problem.c itemset.c queue.c aheap.c sgraph.c stdlib2.c vec.c
	gcc -O3 -Os -s -o $@  mace.c

clean:
	rm -f $(OBJ)

.PHONY:
	MACE clean
