#CC=gfortran
CC=f95
FILE=ltrain.f

#all: ltrain_new32 ltrain_new ltrain_new_debug
all: ltrain_new32 ltrain_new ltrain_new_debug amtrain mtrain

ltrain_new: $(FILE)
	$(CC)  $(FILE) -v -O2 -o ltrain_new 

ltrain_new_debug: $(FILE)
	$(CC)  $(FILE) -v -O2 -g -o ltrain_new_debug

ltrain_new32: $(FILE)
	$(CC) -v $(FILE) -O0 -g -o ltrain_new32 -m32 -Bstatic -std=legacy
#
# the followin requires the source file from Tatiana's public, that isn't a part of the project...
# !! EXLUDED from the all target, does not run unless asked!
ltrain_alt32: ltrain_01.f
	$(CC) -v ltrain_01.f -O0 -o ltrain_alt32


mtrain: mtrain.f
	$(CC) -v mtrain.f -v -O2 -o mtrain

# Current version
amtrain: amtrain.f
	$(CC) -v amtrain.f -v -O2 -o amtrain 




