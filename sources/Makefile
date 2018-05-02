# -----------------------------------------------------------------------------------------------------------
#  variables
# -----------------------------------------------------------------------------------------------------------
# nom projet
PROJECT = sph

#
SOURCES = main.f95 sph.f95 donnees.f95 math.f95

OBJ_DIR = ../objets
RUN_DIR = ../run

OBJECTS = $(wildcard $(OBJ_DIR)/*.o)
MODULES = $(wildcard $(OBJ_DIR)/*.mod)

FC      = gfortran
FLAGS   = -fbounds-check -J $(OBJ_DIR)
EXEC    = truc
OTHER   = ../entrees/constantes
PDF     =

#VPATH = $(OBJ_DIR)



# -----------------------------------------------------------------------------------------------------------
#  compilation rules
# -----------------------------------------------------------------------------------------------------------
# exécutable
$(RUN_DIR)/$(EXEC) : $(OBJ_DIR)/main.o \
	                 $(OBJ_DIR)/math.o \
					 $(OBJ_DIR)/donnees.o \
					 $(OBJ_DIR)/sph.o
	$(FC) $(FLAGS) $(OBJ_DIR)/*.o -o $(RUN_DIR)/$(EXEC)

# programme principal
$(OBJ_DIR)/main.o : main.f95 \
	                $(OBJ_DIR)/math.mod \
					$(OBJ_DIR)/donnees.mod \
					$(OBJ_DIR)/sph.mod
	$(FC) $(FLAGS) -c main.f95 -o $(OBJ_DIR)/main.o

# modules
$(OBJ_DIR)/math.o $(OBJ_DIR)/math.mod : math.f95
	$(FC) $(FLAGS) -c math.f95 -o $(OBJ_DIR)/math.o

$(OBJ_DIR)/donnees.o $(OBJ_DIR)/donnees.mod : donnees.f95 $(OBJ_DIR)/math.mod
	$(FC) $(FLAGS) -c donnees.f95 -o $(OBJ_DIR)/donnees.o

$(OBJ_DIR)/sph.o $(OBJ_DIR)/sph.mod : sph.f95 $(OBJ_DIR)/math.mod
	$(FC) $(FLAGS) -c sph.f95 -o $(OBJ_DIR)/sph.o



# -----------------------------------------------------------------------------------------------------------
# phony targets
# -----------------------------------------------------------------------------------------------------------
# exécuter l'exécutable (utile pour utiliser F5 dans vim)
.PHONY : make
make :
	@ cd $(RUN_DIR); ./$(EXEC)

# supprimer les fichiers objet et l'exécutable s'ils existent
.PHONY : clean
clean :
	rm -f $(OBJ_DIR)/* $(RUN_DIR)/*

# effacer le contenu des dossiers d'entrees et de sorties
.PHONY : del
del :
	rm -f ../sorties/* ../graphes/*

# ouvrir les fichiers du projet dans des onglets de vim
.PHONY : open
open :
	@ vim -p $(SOURCES) Makefile Plot.gnu $(OTHER)

# tout compiler et lancer gdb (segmentation fault)
.PHONY : gdb
gdb :
	$(FC) -g $(SOURCES) -o $(EXEC) && gdb ./$(EXEC)

# clean et tarer le dossier
.PHONY : tar
tar :
	make clean
	cd ../..; tar -zcvf $(PROJECT).tar.gz $(PROJECT)

# sauvegarder ancienne version
.PHONY : save
save :
	make clean
	cd ../..; cp -r $(PROJECT) old_$(PROJECT)

.PHONY : pdf
pdf :
	@ xdg-open $(PDF)

#
.PHONY : clean
coffe :
	@ echo "  (\n   )\n c[]"
