PLANTUML=~/bin/plantuml.jar

UMLFILES= \
InputAdapterClass.txt \
ObjectiveClass.txt \
ProblemClass.txt \
IdentifierMapClass.txt

PNGFILES=$(UMLFILES:.txt=.png)

all: $(PNGFILES)

%.png: %.txt
	java -jar $(PLANTUML) $<
