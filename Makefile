CC=gcc
CFLAGS=-I.
DEPS = utils.h quad_model.h sensors_model.h quad_parameters.h px4_sim_communication.h px4_quad_sim.h
LIBS=-lm

_OBJ = sensors_model.o utils.o quad_model.o px4_sim_communication.o px4_quad_sim.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

ODIR=obj
DIR_GUARD=@mkdir -p $(@D)

$(ODIR)/%.o: %.c $(DEPS)
	$(DIR_GUARD)
	$(CC) -c -o $@ $< $(CFLAGS)

px4_quad_sim: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -rf $(ODIR)
