default:
	@echo "`make scs`?"

scs:
	unoconv -f pdf StepControlStrategy.pptx
	convert StepControlStrategy.pdf StepControlStrategy.png

isd:
	unoconv -f pdf IntegratorStepperDesign.pptx
	convert IntegratorStepperDesign.pdf IntegratorStepperDesign.png

bdf:
	unoconv -f pdf BDF.pptx
	convert BDF.pdf BDF.png

ib:
	unoconv -f pdf InterpolationBuffer.pptx
	convert InterpolationBuffer.pdf InterpolationBuffer.png

clean:
	rm -f BDF.pdf BDF.png \
	  StepControlStrategy.pdf StepControlStrategy.png \
	  IntegratorStepperDesign.pdf IntegratorStepperDesign.png \
	  InterpolationBuffer.pdf InterpolationBuffer.png
