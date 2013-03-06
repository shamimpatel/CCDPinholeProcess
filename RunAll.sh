./CCDPinholeProcessMain || exit 1
sleep 1s
./PlotCCD.sh
python ProduceSpectrum.py DiffractResultsPostPinhole.txt || exit 1
sleep 1s
./PlotSpectrum.sh
sleep 1s
./CorrelationPlot.sh