TARGET=quadratic

AUX=../ML_learning.bib Lagrange.eps

$(TARGET).pdf: $(TARGET).tex $(AUX)
	latex $(TARGET)
	bibtex $(TARGET)
	latex $(TARGET)
	latex $(TARGET)
	dvipdf $(TARGET)
