./CCDPinholeProcessMain || exit 1
sleep 0.5s
./PlotCCD.sh
#python ProduceSpectrum.py DiffractResultsPostPinhole.txt || exit 1
python ProduceBinnedSpectrum.py
sleep 0.5s
./PlotSpectrum.sh
sleep 0.5s
./CorrelationPlot.sh