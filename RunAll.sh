./CCDPinholeProcessMain || exit 1
./PlotCCD.sh
python ProduceSpectrum.py DiffractResultsPostPinhole.txt || exit 1
./PlotSpectrum.sh
./CorrelationPlot.sh