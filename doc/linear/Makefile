TARGET=linear

AUX=../ML_learning.bib

$(TARGET).pdf: $(TARGET).tex $(AUX)
	latex $(TARGET)
	#bibtex $(TARGET)
	#latex $(TARGET)
	latex $(TARGET)
	dvipdf $(TARGET)
