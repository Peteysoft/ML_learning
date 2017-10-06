TARGET=Bayesian

AUX=logistic1.eps roc.eps sample_class.eps

$(TARGET).pdf: $(TARGET).tex $(AUX)
	latex $(TARGET)
	latex $(TARGET)
	dvipdf $(TARGET)
