OBJ = dafor.o

dafor: $(OBJ) dafor.mk
	g77 -o dafor $(OBJ)
