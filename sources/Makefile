# ===========================================================================================================
# Variables.
# ===========================================================================================================
# Project name.
PROJECT = sources

# Structure.
OBJ_DIR = ../objects
RUN_DIR = ../run
EXEC    = truc

# Files to edit.
SOURCES = main.f95 donnees.f95 sph.f95 tsup.f95 kernel.f95 var.f95 math.f95
MODULES = donnees.mod sph.mod tsup.mod kernel.mod var.mod math.mod
OBJECTS = main.o donnees.o sph.o tsup.o kernel.o var.o math.o
OTHER   = Makefile ../entrees/constantes Plot.gnu
PDF     =

# Compilator informations.
FC      = gfortran
FLAGS   = -fbounds-check -J $(OBJ_DIR) -cpp



# ===========================================================================================================
# Compilationrules.
# ===========================================================================================================
# Exécutable
$(RUN_DIR)/$(EXEC) : $(OBJ_DIR)/main.o \
                     $(OBJ_DIR)/donnees.o \
                     $(OBJ_DIR)/sph.o \
                     $(OBJ_DIR)/tsup.o \
                     $(OBJ_DIR)/kernel.o \
                     $(OBJ_DIR)/var.o \
                     $(OBJ_DIR)/math.o
	$(FC) $(FLAGS) $(OBJ_DIR)/*.o -o $(RUN_DIR)/$(EXEC)


# Programme principal
$(OBJ_DIR)/main.o : main.f95 \
                    $(OBJ_DIR)/sph.mod \
                    $(OBJ_DIR)/donnees.mod
	$(FC) $(FLAGS) -c main.f95 -o $(OBJ_DIR)/main.o


# Modules
$(OBJ_DIR)/donnees.o $(OBJ_DIR)/donnees.mod : donnees.f95 \
                                              $(OBJ_DIR)/sph.mod
	$(FC) $(FLAGS) -c donnees.f95 -o $(OBJ_DIR)/donnees.o

$(OBJ_DIR)/sph.o $(OBJ_DIR)/sph.mod : sph.f95 \
                                      $(OBJ_DIR)/tsup.mod
	$(FC) $(FLAGS) -c sph.f95 -o $(OBJ_DIR)/sph.o

$(OBJ_DIR)/tsup.o $(OBJ_DIR)/tsup.mod : tsup.f95 \
                                        $(OBJ_DIR)/var.mod \
                                        $(OBJ_DIR)/kernel.mod
	$(FC) $(FLAGS) -c tsup.f95 -o $(OBJ_DIR)/tsup.o

$(OBJ_DIR)/kernel.o $(OBJ_DIR)/kernel.mod : kernel.f95 \
                                            $(OBJ_DIR)/var.mod
	$(FC) $(FLAGS) -c kernel.f95 -o $(OBJ_DIR)/kernel.o

$(OBJ_DIR)/var.o $(OBJ_DIR)/var.mod : var.f95 \
                                      $(OBJ_DIR)/math.mod
	$(FC) $(FLAGS) -c var.f95 -o $(OBJ_DIR)/var.o

$(OBJ_DIR)/math.o $(OBJ_DIR)/math.mod : math.f95
	$(FC) $(FLAGS) -c math.f95 -o $(OBJ_DIR)/math.o



# ===========================================================================================================
# phony targets
# ===========================================================================================================
# exécuter l'exécutable (utile pour utiliser F5 dans vim)
.PHONY : make
make :
	@ cd $(RUN_DIR); ./$(EXEC)

# supprimer les fichiers objet et l'exécutable s'ils existent
.PHONY : clean
clean :
	rm -f $(OBJ_DIR)/*.mod $(OBJ_DIR)/*.o $(RUN_DIR)/$(EXEC)

# effacer le contenu des dossiers d'entrees et de sorties
.PHONY : del
del :
	rm -f sorties/*.dat graphes/*.eps graphes/*.ps graphes/*.png

# ouvrir les fichiers du projet dans des onglets de vim
.PHONY : open
open :
	@ vim $(SOURCES) $(OTHER)

# tout compiler et lancer gdb (segmentation fault)
.PHONY : gdb
gdb :
	$(FC) -g $(SOURCES) -o $(EXEC) && (cd $(RUN_DIR); gdb ./$(EXEC))

# clean et tarer le dossier
.PHONY : tar
tar :
	make clean
	cd ..; tar -zcvf $(PROJECT).tar.gz $(PROJECT)

# sauvegarder ancienne version
.PHONY : save
save :
	make clean
	cd ..; cp -r $(PROJECT) old_$(PROJECT)

.PHONY : pdf
pdf :
	@ xdg-open $(PDF)

#
.PHONY : clean
coffe :
	@ echo "  (\n   )\n c[]"
