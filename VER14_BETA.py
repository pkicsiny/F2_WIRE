# Logistical Shit
import sys
sys.path.append('/usr/local/facet/tools/python')
sys.path.append('/usr/local/facet/tools/python/F2_pytools')
from F2_pytools.mags import set_magnets
from pathlib import Path

# SLAC Shit
from physicselog import submit_entry
from epics import get_pv, caget, caput

# Time Shit
from datetime import datetime
from time import sleep, monotonic

# imports for the QT Shit
from PyQt5 import QtWidgets, uic
from PyQt5.QtWidgets import QVBoxLayout
from PyQt5.QtGui import QDoubleValidator, QIntValidator

# Math/Numerical Shit
import numpy as np
import pandas as pd
from numpy.ma.core import shape
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, peak_widths, savgol_filter
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import elog_client as ELOG
########################################################################################################################

#ELOG.post(
#title = title
#body=""
#log=""
#logbook
#file_paths)


# SCRIPT_DIR = Path(__file__).resolve().parent
SCRIPT_DIR = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
# UI_FILENAME = "WireScanner_FirstPASS_6_gpt_light.ui"
UI_FILENAME = "VER14_NWA_BETA.ui"
# "WireScanner_FirstPASS_6_gpt_light.ui"
UI_FALLBACK = Path(
    "/Users/mparker/Library/CloudStorage/OneDrive-SLACNationalAcceleratorLaboratory/ACTIVE WIRE_SCANNER FILE/" + UI_FILENAME)


# UI_FALLBACK = Path("pybox/ProjBox/Rev1_Working Folder")
# UI_FALLBACK = Path("/full/path/to/Rev1_Working Folder/WireScanner_FirstPASS_6_gpt_light.ui")


def find_ui_file():
    """Locate the Qt Designer file without depending on the launch directory."""
    local_ui = SCRIPT_DIR / UI_FILENAME
    if local_ui.is_file():
        return local_ui
    if UI_FALLBACK.is_file():
        return UI_FALLBACK
    raise FileNotFoundError(f"Could not find {UI_FILENAME!r} beside this script or at " f"{UI_FALLBACK}")


########################################################################################################################


class MyGui(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        uic.loadUi(str(find_ui_file()), self)

        self.tabWidget.removeTab(self.tabWidget.indexOf(self.EmittanceMatching))
        #elf.tabWidget.removeTab(self.tabWidget.indexOf(self.AllWirePlots))
        #self.tabWidget.removeTab(self.tabWidget.indexOf(self.L3EnsemblePlots))

        # Quad scan tab: only allow numbers in the boxes
        self.QuadKMinBox.setValidator(QDoubleValidator(-1000.0, 1000.0, 4))
        self.QuadKMaxBox.setValidator(QDoubleValidator(-1000.0, 1000.0, 4))
        self.QuadStepsBox.setValidator(QIntValidator(1, 1000))

        # Quad scan tab: build the k array whenever a value is entered
        self.quad_k_array = None
        self.QuadKMinBox.editingFinished.connect(self.update_quad_k_array)
        self.QuadKMaxBox.editingFinished.connect(self.update_quad_k_array)
        self.QuadStepsBox.editingFinished.connect(self.update_quad_k_array)
        self.quad_scan_fitted_rms_dict = {}
        self.quad_scan_fitted_rms_error_dict = {}

        self.abort_requested = False
        self.axis = None
        self.last_scan = None
        self.multi_scan_results = {
            ("L2", "X"): {},
            ("L2", "Y"): {},
            ("L3", "X"): {},
            ("L3", "Y"): {}, }

        for progress_bar in [self.progressBar_1, self.progressBar_2, self.progressBar_3, self.progressBar_4, self.progressBar_5, self.progressBar_6,]:
            progress_bar.setValue(0)



        # *********
        def make_plot_frame(frame):
            figure = Figure()
            canvas = FigureCanvas(figure)
            axes = figure.add_subplot(111)

            layout = QVBoxLayout(frame)
            layout.setContentsMargins(0, 0, 0, 0)
            layout.setSpacing(0)
            layout.addWidget(canvas)

            return figure, canvas, axes

        def make_1x4_plot_frame(frame):
            figure = Figure()
            canvas = FigureCanvas(figure)
            axes = figure.subplots(1, 4)

            layout = frame.layout()
            if layout is None:
                layout = QVBoxLayout(frame)
                layout.setContentsMargins(0, 0, 0, 0)
                layout.setSpacing(0)

            layout.addWidget(canvas)
            return figure, canvas, axes
        ## Plotting Frames:
        ### L2 X Frames
        self.frame_7_figure, self.frame_7_canvas, self.frame_7_axes = make_plot_frame(self.frame_7)
        ### L2 Y Frames
        self.frame_13_figure, self.frame_13_canvas, self.frame_13_axes = make_plot_frame(self.frame_13)
        ### L3 X Frames
        self.l3XPlot1_figure, self.l3XPlot1_canvas, self.l3XPlot1_axes = make_plot_frame(self.l3XPlot1)
        ### L3 Y Frames
        self.l3YPlot1_figure, self.l3YPlot1_canvas, self.l3YPlot1_axes = make_plot_frame(self.l3YPlot1)

        ### Single wire & IP X|Y Frames
        self.SingleWirePlot_figure, self.SingleWirePlot_canvas, self.SingleWirePlot_axes = make_plot_frame(self.SingleWirePlot_X)
        self.single_scan_results = {}
        #self.Xaxis_Label.hide()
        #self.Yaxis_Label.hide()

        # quad scan frame
        self.QuadScanPlot_figure, self.QuadScanPlot_canvas, self.QuadScanPlot_axes = make_plot_frame(
        self.QuadScanPlot)

        ################################################################################################################
        # Quad Scan Wire scan Tab
        ################################################################################################################

        # Scan Buttons:
        self.QuadScanXYpush.clicked.connect(lambda: self.quad_scan_all(axis="XY"))
        self.QuadScanXpush.clicked.connect( lambda: self.quad_scan_all(axis= "X"))
        self.QuadScanYpush.clicked.connect( lambda: self.quad_scan_all(axis= "Y"))

        ################################################################################################################
        # Single Wire scan Tab
        ################################################################################################################
        # Scan Buttons:
        self.ScanXYpush.clicked.connect(lambda: self.scan_all(area="IP", axis="XY"))
        self.ScanXpush.clicked.connect(lambda: self.singlesscan(axis="X", progressBar=self.progressBar_1, ))
        self.ScanYpush.clicked.connect(lambda: self.singlesscan(axis="Y", progressBar=self.progressBar_1, ))
        self.AbortPush_2.clicked.connect(self.abort_scan)
        # Park Buttons:
        self.ParkPush.clicked.connect(lambda: self.PARKXY())  # Parks to X or Y whatever is selected.
        self.ParkXpush.clicked.connect(lambda: self.parkFunc(axis="X"))
        self.ParkYpush.clicked.connect(lambda: self.parkFunc(axis="Y"))
        # Elog Button:
        self.ELOG.clicked.connect(lambda: self.elog_plot())
        # Clear Buttons
        self.ClearX.clicked.connect(lambda: self.clear_single_plot("X"))
        self.ClearY.clicked.connect(lambda: self.clear_single_plot("Y"))

        # Wire Dropdown:
        self.WireDropdown.clear()
        self.WireDropdown.addItems(
            [
                # L2
                "WIRE:LI11:444", "WIRE:LI11:614", "WIRE:LI11:744", "WIRE:LI12:214",
                # L3
                "WIRE:LI18:944", "WIRE:LI19:144", "WIRE:LI19:244", "WIRE:LI19:344",
                # IP
                "WIRE:LI20:3179",
                # "WIRE:LI20:3229",
            ])
        self.rms_axis_indices = {"X": 0, "Y": 1}
        self.rms_wire_indices = {self.WireDropdown.itemText(index): index for index in range(self.WireDropdown.count())}
        self.rms_result_names = ("raw_rms", "fitted_rms", "raw_rms_error", "fitted_rms_error")
        self.RMSdata = np.full((2, self.WireDropdown.count(), 4), np.nan)
        self.emittance_results = {}

        self.WireDropdown.setCurrentText("WIRE:LI20:3179")
        # self.WireDropdown.currentTextChanged.connect(lambda: self.load_dummy_DATA_PLOTS())
        self.WireDropdown.currentTextChanged.connect(self.wire_selection_changed)

        ################################################################################################################
        # Multi L2 X Wire scan Tab
        ################################################################################################################
        # Scan Buttons:
        self.l2XScan1Push.clicked.connect(
            lambda: self.scan(sec="LI11", loc="444", axis="X", plot_axes=self.frame_7_axes,
                              plot_canvas=self.frame_7_canvas, progressBar=self.progressBar_2))
        self.l2XScan2Push.clicked.connect(
            lambda: self.scan(sec="LI11", loc="614", axis="X", plot_axes=self.frame_7_axes,
                              plot_canvas=self.frame_7_canvas, progressBar=self.progressBar_2))
        self.l2XScan3Push.clicked.connect(
            lambda: self.scan(sec="LI11", loc="744", axis="X", plot_axes=self.frame_7_axes,
                              plot_canvas=self.frame_7_canvas, progressBar=self.progressBar_2))
        self.l2XScan4Push.clicked.connect(
            lambda: self.scan(sec="LI12", loc="214", axis="X", plot_axes=self.frame_7_axes,
                              plot_canvas=self.frame_7_canvas, progressBar=self.progressBar_2))
        # Logistical Buttons
        self.l2XElogPush.clicked.connect(lambda: self.elog_plot())
        self.l2XAbortPush.clicked.connect(self.abort_scan)
        self.l2XScanAllPush.clicked.connect(lambda: self.scan_all(area="L2", axis="X"))

        # CLEAR Buttons

        ################################################################################################################
        # Multi L2 Y Wire scan Tab
        ################################################################################################################
        self.yScan1Push.clicked.connect(lambda: self.scan(sec="LI11", loc="444", axis="Y", plot_axes=self.frame_13_axes,
                                                          plot_canvas=self.frame_13_canvas,
                                                          progressBar=self.progressBar_3))
        self.yScan2Push.clicked.connect(lambda: self.scan(sec="LI11", loc="614", axis="Y", plot_axes=self.frame_13_axes,
                                                          plot_canvas=self.frame_13_canvas,
                                                          progressBar=self.progressBar_3))
        self.yScan3Push.clicked.connect(lambda: self.scan(sec="LI11", loc="744", axis="Y", plot_axes=self.frame_13_axes,
                                                          plot_canvas=self.frame_13_canvas,
                                                          progressBar=self.progressBar_3))
        self.yScan4Push.clicked.connect(lambda: self.scan(sec="LI12", loc="214", axis="Y", plot_axes=self.frame_13_axes,
                                                          plot_canvas=self.frame_13_canvas,
                                                          progressBar=self.progressBar_3))
        self.yElogPush.clicked.connect(lambda: self.elog_plot())
        self.yScanAllPush.clicked.connect(lambda: self.scan_all(area="L2", axis="Y"))
        self.yAbortPush.clicked.connect(self.abort_scan)

        ################################################################################################################
        # Multi L3 X Wire scan Tab
        ################################################################################################################
        self.l3XScan1Push.clicked.connect(
            lambda: self.scan(sec="LI18", loc="944", axis="X", plot_axes=self.l3XPlot1_axes,
                              plot_canvas=self.l3XPlot1_canvas, progressBar=self.progressBar_4))
        self.l3XScan2Push.clicked.connect(
            lambda: self.scan(sec="LI19", loc="144", axis="X", plot_axes=self.l3XPlot1_axes,
                              plot_canvas=self.l3XPlot1_canvas, progressBar=self.progressBar_4))
        self.l3XScan3Push.clicked.connect(
            lambda: self.scan(sec="LI19", loc="244", axis="X", plot_axes=self.l3XPlot1_axes,
                              plot_canvas=self.l3XPlot1_canvas, progressBar=self.progressBar_4))
        self.l3XScan4Push.clicked.connect(
            lambda: self.scan(sec="LI19", loc="344", axis="X", plot_axes=self.l3XPlot1_axes,
                              plot_canvas=self.l3XPlot1_canvas, progressBar=self.progressBar_4))
        self.l3XElogPush.clicked.connect(lambda: self.elog_plot())
        self.l3XScanAllPush.clicked.connect(lambda: self.scan_all(area="L3", axis="X"))
        self.l3XAbortPush.clicked.connect(self.abort_scan)

        ################################################################################################################
        # Multi L3 Y Wire scan Tab
        ################################################################################################################
        self.l3YScan1Push.clicked.connect(
            lambda: self.scan(sec="LI18", loc="944", axis="Y", plot_axes=self.l3YPlot1_axes,
                              plot_canvas=self.l3YPlot1_canvas, progressBar=self.progressBar_5))
        self.l3YScan2Push.clicked.connect(
            lambda: self.scan(sec="LI19", loc="144", axis="Y", plot_axes=self.l3YPlot1_axes,
                              plot_canvas=self.l3YPlot1_canvas, progressBar=self.progressBar_5))
        self.l3YScan3Push.clicked.connect(
            lambda: self.scan(sec="LI19", loc="244", axis="Y", plot_axes=self.l3YPlot1_axes,
                              plot_canvas=self.l3YPlot1_canvas, progressBar=self.progressBar_5))
        self.l3YScan4Push.clicked.connect(
            lambda: self.scan(sec="LI19", loc="344", axis="Y", plot_axes=self.l3YPlot1_axes,
                              plot_canvas=self.l3YPlot1_canvas, progressBar=self.progressBar_5))
        self.l3YElogPush.clicked.connect(lambda: self.elog_plot())
        self.l3YScanAllPush.clicked.connect(lambda: self.scan_all(area="L3", axis="Y"))
        self.l3YAbortPush.clicked.connect(self.abort_scan)

        ################################################################################################################
        # Emit/Match Tab
        ################################################################################################################
        # self.l2CalcEmitXPush.clicked.connect(lambda checked=False: self.run_emittance_calculation("L2", "X"))
        # self.l2CalcEmitYPush.clicked.connect(lambda checked=False: self.run_emittance_calculation("L2", "Y"))
        # self.l3CalcEmitXPush.clicked.connect(lambda checked=False: self.run_emittance_calculation("L3", "X"))
        # self.l3CalcEmitYPush.clicked.connect(lambda checked=False: self.run_emittance_calculation("L3", "Y"))

    def clear_single_plot(self, axis):
        self.single_scan_results.pop(axis, None)
        self.update_single_subplots()
        print(f"Cleared {axis} plot")

    def abort_scan(self):
        self.abort_requested = True
        print("Scan aborted")

    def get_multi_plot(self, area, axis):
        area = area.upper()
        axis = axis.upper()

        plot_map = {
            ("L2", "X"): (
                self.frame_7_figure,
                self.frame_7_canvas,
            ),
            ("L2", "Y"): (
                self.frame_13_figure,
                self.frame_13_canvas,
            ),
            ("L3", "X"): (
                self.l3XPlot1_figure,
                self.l3XPlot1_canvas,
            ),
            ("L3", "Y"): (
                self.l3YPlot1_figure,
                self.l3YPlot1_canvas,
            ),
        }

        key = (area, axis)

        if key not in plot_map:
            raise ValueError(
                f"Unsupported plot selection: area={area}, axis={axis}"
            )

        return plot_map[key]

    def update_multi_subplots(self, area, axis):
        area = area.upper()
        axis = axis.upper()
        key = (area, axis)

        figure, canvas = self.get_multi_plot(area, axis)
        scan_results = self.multi_scan_results[key]
        number_of_plots = len(scan_results)

        figure.clear()

        if number_of_plots == 0:
            canvas.draw_idle()
            return
        elif number_of_plots == 1:
            rows, columns = 1, 1
        elif number_of_plots == 2:
            rows, columns = 1, 2
        elif number_of_plots == 3:
            rows, columns = 1, 3
        else:
            rows, columns = 2, 2

        axes = figure.subplots(rows, columns, squeeze=False, ).flatten()

        for plot_axis, (wire_name, scan_data) in zip(
                axes,
                scan_results.items(),
        ):
            scan_data["finder"].plot_beam_rms(beam_result=scan_data["result"], error_result=scan_data["error_result"], show=False, ax=plot_axis,)

            plot_axis.set_title(f"{wire_name} {axis}")

        # Remove any unused subplot positions.
        for unused_axis in axes[number_of_plots:]:
            unused_axis.remove()

        figure.tight_layout()
        canvas.draw_idle()

    def update_single_subplots(self):
        figure = self.SingleWirePlot_figure
        figure.clear()

        scanned_axes = [
            axis for axis in ("X", "Y")
            if axis in self.single_scan_results
        ]

        if len(scanned_axes) == 1:
            plot_axes = [figure.add_subplot(1, 1, 1)]
        else:
            plot_axes = figure.subplots(1, 2, squeeze=False).flatten()

        for axis, ax in zip(scanned_axes, plot_axes):
            scan = self.single_scan_results[axis]

            scan["finder"].plot_beam_rms(
                beam_result=scan["result"],
                error_result=scan["error_result"],
                show=False,
                ax=ax,
            )

            ax.set_title(f"{axis}-Axis")

        figure.tight_layout()
        self.SingleWirePlot_canvas.draw_idle()

    def quad_singlesscan(self, axis=None, progressBar=None):

        ####################       
        # init wire object #
        ####################
 
        self.abort_requested = False
        axis = "X" if axis is None else axis.upper()

        self.axis = axis
        self.area = "IP"
        self.sec = "LI20"
        self.loc = "3179"
        self.wirePV = f"WIRE:{self.sec}:{self.loc}"

        wire = wireFunk(wirePV=self.wirePV, area=self.area, axis=axis)

        #####################
        # init progress bar #
        #####################

        if progressBar is not None:
            progressBar.setValue(0)

        def update_scan_progress(value):
            if progressBar is not None:
                progressBar.setValue(value)
                QtWidgets.QApplication.processEvents()
            return self.abort_requested

        ####################
        # do the wire scan #
        ####################

        x, PMTdata = wire.scanfunk(progress_callback=update_scan_progress)

        if self.abort_requested:
            print(f"{self.axis} scan aborted")
            return

        ################
        # process data #
        ################

        smooth_signal = PMTdata.copy()
        finder = PF_FINDER(x, smooth_signal, self.wirePV)
        finder.find_peaks_and_widths()

        try:
            result = finder.calculate_beam_rms(raw_y=PMTdata, smoothed_y=smooth_signal, )
            error_result = finder.calculate_single_scan_error(raw_y=PMTdata, smoothed_y=smooth_signal, n_bootstrap=3, )
        except RuntimeError:
            print("WARNING: Gaussian fit failed; using raw beam values")

            result = finder.calculate_beam_rms(raw_y=PMTdata, smoothed_y=smooth_signal, fit_asym_gauss=False, )

            result["fitted_rms"] = result["raw_rms"]
            result["corrected_fitted_rms"] = result["corrected_raw_rms"]
            result["fitted_centroid"] = result["beam_centroid"]
            result["fitted_skew"] = result["beam_skew"]
            result["fitted_kurtosis"] = result["beam_kurtosis"]

            error_result = {"raw_rms_error": 0.0, "fitted_rms_error": 0.0, }

        self.single_scan_results[axis] = {"finder": finder, "result": result, "error_result": error_result, }

        ###############
        # plot on gui #
        ###############

        self.update_single_subplots()
        self.SingleWirePlot_canvas.draw_idle()

        # Save the figure from this scan.
        # Save the current single-wire figure.
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
        self.outfile = Path(f"WIRE_{self.sec}_{self.loc}_{axis}_{timestamp}.png")

        self.SingleWirePlot_figure.savefig(self.outfile, dpi=150, bbox_inches="tight", )

        # Remember the scan for the Elog button.
        self.last_scan = {
            "wire": f"WIRE:{self.sec}:{self.loc}",
            "axis": axis,
            "figure": self.SingleWirePlot_figure,
            "outfile": self.outfile,
            "raw_rms": result["raw_rms"],
            "raw_rms_error": error_result["raw_rms_error"],
            "fitted_rms": result["fitted_rms"],
            "fitted_rms_error": error_result["fitted_rms_error"],
            "raw_skew": result["beam_skew"],
            "raw_kurtosis": result["beam_kurtosis"],
            "fitted_skew": result["fitted_skew"],
            "fitted_kurtosis": result["fitted_kurtosis"],
        }

        print(f"Figure saved: {self.outfile.resolve()}")
        self.elog_plot()
        print("Sent to Facet Elog")

        if self.axis == "X":
            self.InnerposX.setText(str(wire.OGlower_pv.get()))
            self.OutterposX.setText(str(wire.OGupper_pv.get()))
            self.X_Size.setText(f'{result["fitted_rms"]:.3f} ± {error_result["fitted_rms_error"]:.3f}')
            # if self.sec == "LI20" and self.loc == "3179":
            # ep.caput("WIRE:LI20:3179:XRMS", result["raw_rms"], )

        elif self.axis == "Y":
            self.InnerposY.setText(str(wire.OGlower_pv.get()))
            self.OutterposY.setText(str(wire.OGupper_pv.get()))
            self.Y_Size.setText(f'{result["fitted_rms"]:.3f} ± {error_result["fitted_rms_error"]:.3f}')


        ##################
        # append to dict #
        ##################

        if axis not in self.quad_scan_fitted_rms_dict.keys():
            self.quad_scan_fitted_rms_dict[axis] = [self.last_scan["fitted_rms"]]
            self.quad_scan_fitted_rms_error_dict[axis] = [self.last_scan["fitted_rms_error"]]
        else:
            self.quad_scan_fitted_rms_dict[axis].append(self.last_scan["fitted_rms"])
            self.quad_scan_fitted_rms_error_dict[axis].append(self.last_scan["fitted_rms_error"])


    def singlesscan(self, sec=None, loc=None, axis=None, progressBar=None):
        self.abort_requested = False
        axis = "X" if axis is None else axis.upper()
        self.axis = axis

        if sec is None or loc is None:
            _, selected_sec, selected_loc = (self.WireDropdown.currentText().split(":"))
            sec = selected_sec if sec is None else sec
            loc = selected_loc if loc is None else loc

        if sec == "LI11" or sec == "LI12":
            self.area = "L2"
        elif sec == "LI18" or sec == "LI19":
            self.area = "L3"
        elif sec == "LI20":
            self.area = "IP"

        self.sec = sec
        self.loc = loc
        self.wirePV = f"WIRE:{self.sec}:{self.loc}"
        wire = wireFunk(wirePV=self.wirePV, area=self.area, axis=axis)

        if progressBar is not None:
            progressBar.setValue(0)

        def update_scan_progress(value):
            if progressBar is not None:
                progressBar.setValue(value)
                QtWidgets.QApplication.processEvents()
            return self.abort_requested

        x, PMTdata = wire.scanfunk(progress_callback=update_scan_progress)

        if self.abort_requested:
            print(f"{self.axis} scan aborted")
            return

        smooth_signal = PMTdata.copy()
        # finder = PF_FINDER(x,wire=self.wire, y=smooth_signal)
        finder = PF_FINDER(x, smooth_signal, self.wirePV)
        # finder.wirePV = self.wirePV
        finder.find_peaks_and_widths()

        # result = finder.calculate_beam_rms(raw_y=PMTdata, smoothed_y=smooth_signal, )
        # error_result = finder.calculate_single_scan_error(raw_y=PMTdata, smoothed_y=smooth_signal, n_bootstrap=100, )
        # finder.plot_beam_rms(beam_result=result, error_result=error_result, show=False, ax=plot_axes, )
        try:
            result = finder.calculate_beam_rms(raw_y=PMTdata, smoothed_y=smooth_signal, )
            error_result = finder.calculate_single_scan_error(raw_y=PMTdata, smoothed_y=smooth_signal, n_bootstrap=3, )
        except RuntimeError:
            print("WARNING: Gaussian fit failed; using raw beam values")

            result = finder.calculate_beam_rms(raw_y=PMTdata, smoothed_y=smooth_signal, fit_asym_gauss=False, )

            result["fitted_rms"] = result["raw_rms"]
            result["corrected_fitted_rms"] = result["corrected_raw_rms"]
            result["fitted_centroid"] = result["beam_centroid"]
            result["fitted_skew"] = result["beam_skew"]
            result["fitted_kurtosis"] = result["beam_kurtosis"]

            error_result = {"raw_rms_error": 0.0, "fitted_rms_error": 0.0, }

        # finder.plot_beam_rms(beam_result=result, error_result=error_result, show=False, ax=plot_axes,)
        self.single_scan_results[axis] = {"finder": finder, "result": result, "error_result": error_result, }

        self.update_single_subplots()
        self.SingleWirePlot_canvas.draw_idle()

        # self.show_single_plot(axis)
        # finder.plot_beam_rms(beam_result=result, error_result=error_result, show=False, ax=plot_axes,)
        # plot_canvas.draw_idle()

        # Save the figure from this scan.
        # Save the current single-wire figure.
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
        self.outfile = Path(f"WIRE_{sec}_{loc}_{axis}_{timestamp}.png")

        self.SingleWirePlot_figure.savefig(self.outfile, dpi=150, bbox_inches="tight", )

        # Remember the scan for the Elog button.
        self.last_scan = {
            "wire": f"WIRE:{sec}:{loc}",
            "axis": axis,
            "figure": self.SingleWirePlot_figure,
            "outfile": self.outfile,
            "raw_rms": result["raw_rms"],
            "raw_rms_error": error_result["raw_rms_error"],
            "fitted_rms": result["fitted_rms"],
            "fitted_rms_error": error_result["fitted_rms_error"],
            "raw_skew": result["beam_skew"],
            "raw_kurtosis": result["beam_kurtosis"],
            "fitted_skew": result["fitted_skew"],
            "fitted_kurtosis": result["fitted_kurtosis"],
        }

        print(f"Figure saved: {self.outfile.resolve()}")
        # self.elog_plot()
        # print("Sent to Facet Elog")
        if self.axis == "X":
            self.InnerposX.setText(str(wire.OGlower_pv.get()))
            self.OutterposX.setText(str(wire.OGupper_pv.get()))
            self.X_Size.setText(f'{result["fitted_rms"]:.3f} ± {error_result["fitted_rms_error"]:.3f}')
            # if self.sec == "LI20" and self.loc == "3179":
            # ep.caput("WIRE:LI20:3179:XRMS", result["raw_rms"], )

        elif self.axis == "Y":
            self.InnerposY.setText(str(wire.OGlower_pv.get()))
            self.OutterposY.setText(str(wire.OGupper_pv.get()))
            self.Y_Size.setText(f'{result["fitted_rms"]:.3f} ± {error_result["fitted_rms_error"]:.3f}')
        # if self.sec == "LI20" and self.loc == "3179":
        # ep.caput("WIRE:LI20:3179:YRMS", result["raw_rms"], )

    def scan(self, sec=None, loc=None, axis=None, plot_axes=None, plot_canvas=None, progressBar=None):
        self.abort_requested = False
        axis = "X" if axis is None else axis.upper()
        self.axis = axis

        if sec is None or loc is None:
            _, selected_sec, selected_loc = self.WireDropdown.currentText().split(":")
            sec = selected_sec if sec is None else sec
            loc = selected_loc if loc is None else loc

        if sec in ("LI11", "LI12"):
            self.area = "L2"
        elif sec in ("LI18", "LI19"):
            self.area = "L3"
        elif sec == "LI20":
            self.area = "IP"
        else:
            raise ValueError(f"Unknown wire section: {sec}")

        self.sec = sec
        self.loc = loc
        self.wirePV = f"WIRE:{self.sec}:{self.loc}"

        wire = wireFunk(wirePV=self.wirePV, area=self.area, axis=axis)

        if progressBar is not None:
            progressBar.setValue(0)

        def update_scan_progress(value):
            if progressBar is not None:
                progressBar.setValue(value)
                QtWidgets.QApplication.processEvents()
            return self.abort_requested

        x, PMTdata = wire.scanfunk(progress_callback=update_scan_progress)

        if self.abort_requested:
            return

        # x, PMTdata = wire.scanfunk(progress_callback=update_scan_progress)
        smooth_signal = savgol_filter(PMTdata, 11, 3) if len(PMTdata) >= 11 else PMTdata.copy()

        finder = PF_FINDER(x, smooth_signal, self.wirePV)
        # finder.wirePV = self.wirePV
        finder.find_peaks_and_widths()

        result = finder.calculate_beam_rms(raw_y=PMTdata, smoothed_y=smooth_signal)
        error_result = finder.calculate_single_scan_error(raw_y=PMTdata, smoothed_y=smooth_signal, n_bootstrap=100)

        wire_name = f"WIRE:{sec}:{loc}"

        # Multi-wire L2/L3 subplot.
        if self.area in ("L2", "L3"):
            key = (self.area, self.axis)
            self.multi_scan_results[key][wire_name] = {"finder": finder, "result": result, "error_result": error_result}
            self.update_multi_subplots(area=self.area, axis=self.axis)
            active_figure, _ = self.get_multi_plot(self.area, self.axis)

        # Single-wire/IP plot.
        else:
            finder.plot_beam_rms(beam_result=result, error_result=error_result, show=False, ax=plot_axes)
            plot_canvas.draw_idle()
            active_figure = plot_axes.figure

        # Store results for the emittance calculation.
        if wire_name in self.rms_wire_indices:
            axis_index = self.rms_axis_indices[axis]
            wire_index = self.rms_wire_indices[wire_name]

            self.RMSdata[axis_index, wire_index, :] = [
                result["raw_rms"],
                result["fitted_rms"],
                error_result["raw_rms_error"],
                error_result["fitted_rms_error"],
            ]

        # Save the active figure.
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
        self.outfile = Path(f"WIRE_{sec}_{loc}_{axis}_{timestamp}.png")
        active_figure.savefig(self.outfile, dpi=150, bbox_inches="tight")

        # Remember the most recent scan.
        self.last_scan = {
            "wire": wire_name,
            "axis": axis,
            "figure": active_figure,
            "outfile": self.outfile,
            "raw_rms": result["raw_rms"],
            "fitted_rms": result["fitted_rms"],
            "raw_rms_error": error_result["raw_rms_error"],
            "fitted_rms_error": error_result["fitted_rms_error"],
            "raw_skew": result["beam_skew"],
            "raw_kurtosis": result["beam_kurtosis"],
            "fitted_skew": result["fitted_skew"],
            "fitted_kurtosis": result["fitted_kurtosis"],
        }

        print(f"Figure saved: {self.outfile.resolve()}")
        # self.elog_plot()
        # print("Sent to Facet Elog")

        if self.sec == "LI20" and self.loc == "3179":
            if self.axis == "X":
                ep.caput("WIRE:LI20:3179:XRMS", result["raw_rms"])
            elif self.axis == "Y":
                ep.caput("WIRE:LI20:3179:YRMS", result["raw_rms"])

    def quad_scan_all(self, axis="XY"):

        """
        only works for area=IP and axis="XY"
        """

        # caputre OG quad value
        quad_pv_root = self.QuadPVDropdown.currentText()
        quad_pv_OG_str = f"{quad_pv_root}:BDES"
        quad_pv_OG = get_pv(quad_pv_OG_str)
        ki_OG = quad_pv_OG.get()

        if self.quad_k_array is None:
            raise ValueError(f"self.quad_k_array is none: {self.quad_k_array}. Set quad scan points first!")
        else:
            
            # new scan: forget the results of previous scans
            self.quad_scan_fitted_rms_dict = {}
            self.quad_scan_fitted_rms_error_dict = {}

            ######################
            # loop over k values #
            ######################

            for ki in self.quad_k_array:
                
                print(f"[quad_scan_all] Initiating quad scan with k1={ki} [kG]")

                #####################
                # set quad strength #
                #####################

                #quad_pv_set_str = f"{quad_pv_root}:BDES"
                quad_pv_get_str = f"{quad_pv_root}:BACT"
                #quad_pv_set = get_pv(quad_pv_set_str)
                quad_pv_get = get_pv(quad_pv_get_str)

                # set quad strength
                set_magnets([quad_pv_root], [float(ki)])
                #quad_pv_set.put(float(ki), wait=True) # quad strength in [kG] 

                # check if properly set
                tol = 1e-1
                ki_readback = quad_pv_get.get()
                while abs(ki_readback-ki) > tol:
                    ki_readback = quad_pv_get.get()
                print(f"[quad_scan_all] Current quad value set to: {ki_readback} [kG]")                   

                #####################
                # perform wire scan #
                #####################

                if axis=="XY":
                    self.quad_singlesscan(axis="X", progressBar=self.progressBar_6)
                    if not self.abort_requested:
                        self.quad_singlesscan(axis="Y", progressBar=self.progressBar_6)
                elif axis=="X":
                    self.quad_singlesscan(axis="X", progressBar=self.progressBar_6)
                elif axis=="Y":
                    self.quad_singlesscan(axis="Y", progressBar=self.progressBar_6)
                else:
                    raise ValueError(f"Invalid value for axis: {axis}. Has to be X, Y or XY!")

                if self.abort_requested:
                    print("[quad_scan_all] aborted, restoring quad")
                    break

                print("[quad_scan_all] wire scan complete")


            # set back OG quad value
            set_magnets([quad_pv_root], [float(ki_OG)])

            # check if properly set
            tol = 1e-1
            ki_OG_readback = quad_pv_get.get()
            while abs(ki_OG_readback-ki_OG) > tol:
                ki_OG_readback = quad_pv_get.get()
            print(f"[quad_scan_all] Quad k value set to: {ki_OG_readback} [kG]")                   

            #######################
            # make plot of k1-rms #
            #######################
                
            self.quad_scan_plot()



    def quad_scan_plot(self):
        """Plot fitted rms vs k1 on the GUI canvas (self.QuadScanPlot_canvas).
        One subplot for an X or a Y scan, two subplots (X | Y) for an XY scan."""
        if self.quad_k_array is None:
            raise ValueError(f"self.quad_k_array is none: {self.quad_k_array}. Set quad scan points first!")

        assert self.quad_scan_fitted_rms_dict.keys() == self.quad_scan_fitted_rms_error_dict.keys(), "fitted rms and fitted rms error dict keys do not agree!"
        axes_keys = [k for k in ("X", "Y") if k in self.quad_scan_fitted_rms_dict]  # ['X', 'Y'] or ['X'] or ['Y']

        # use the figure that lives inside the GUI frame and start from a blank page
        fig = self.QuadScanPlot_figure
        fig.clear()

        if len(axes_keys) == 0:  # nothing measured yet (e.g. aborted during the first point)
            self.QuadScanPlot_canvas.draw_idle()
            return None

        # squeeze=False always returns a 2D array, so ax[i] works for 1 and for 2 subplots
        ax = fig.subplots(1, len(axes_keys), squeeze=False)[0]

        ######################
        # loop over X Y axes #
        ######################

        for i, ax_key in enumerate(axes_keys):

            print(f"[quad_scan_plot] plotting {ax_key} scan data...")
            rms_list = self.quad_scan_fitted_rms_dict[ax_key]
            rms_error_list = self.quad_scan_fitted_rms_error_dict[ax_key]
            assert len(rms_list) == len(rms_error_list), "rms and rms error lists have different length!"

            ########
            # plot #
            ########

            # n < len(k_array) if the scan was aborted
            n = len(rms_list)
            ax[i].errorbar(self.quad_k_array[:n], rms_list, yerr=rms_error_list, c="b", capsize=3)

            #################
            # plot settings #
            #################

            ax[i].set_xlabel("k1 [kG]")
            ax[i].set_ylabel("fitted rms [um]")

            quad_pv_root = self.QuadPVDropdown.currentText()
            ax[i].set_title(f"{quad_pv_root} quad scan, {ax_key} axis wire scan")

        fig.tight_layout()
        self.QuadScanPlot_canvas.draw_idle()

        return fig, ax

    def scan_all(self, area="L2", axis="X"):
        if area == "L2" and axis == "X":
            self.scan(sec="LI11", loc="444", axis="X", plot_axes=self.frame_7_axes, plot_canvas=self.frame_7_canvas,
                      progressBar=self.progressBar_2),
            self.scan(sec="LI11", loc="614", axis="X", plot_axes=self.frame_7_axes, plot_canvas=self.frame_7_canvas,
                      progressBar=self.progressBar_2),
            self.scan(sec="LI11", loc="744", axis="X", plot_axes=self.frame_7_axes, plot_canvas=self.frame_7_canvas,
                      progressBar=self.progressBar_2),
            self.scan(sec="LI12", loc="214", axis="X", plot_axes=self.frame_7_axes, plot_canvas=self.frame_7_canvas,
                      progressBar=self.progressBar_2)
        elif area == "L2" and axis == "Y":
            self.scan(sec="LI11", loc="444", axis="Y", plot_axes=self.frame_13_axes, plot_canvas=self.frame_13_canvas)
            self.scan(sec="LI11", loc="614", axis="Y", plot_axes=self.frame_13_axes, plot_canvas=self.frame_13_canvas)
            self.scan(sec="LI11", loc="744", axis="Y", plot_axes=self.frame_13_axes, plot_canvas=self.frame_13_canvas)
            self.scan(sec="LI12", loc="214", axis="Y", plot_axes=self.frame_13_axes, plot_canvas=self.frame_13_canvas)
        elif area == "L3" and axis == "X":
            self.scan(sec="LI18", loc="944", axis="X", plot_axes=self.l3XPlot1_axes, plot_canvas=self.l3XPlot1_canvas)
            self.scan(sec="LI19", loc="144", axis="X", plot_axes=self.l3XPlot1_axes, plot_canvas=self.l3XPlot1_canvas)
            self.scan(sec="LI19", loc="244", axis="X", plot_axes=self.l3XPlot1_axes, plot_canvas=self.l3XPlot1_canvas)
            self.scan(sec="LI19", loc="344", axis="X", plot_axes=self.l3XPlot1_axes, plot_canvas=self.l3XPlot1_canvas)
        elif area == "L3" and axis == "Y":
            self.scan(sec="LI18", loc="944", axis="Y", plot_axes=self.l3YPlot1_axes, plot_canvas=self.l3YPlot1_canvas)
            self.scan(sec="LI19", loc="144", axis="Y", plot_axes=self.l3YPlot1_axes, plot_canvas=self.l3YPlot1_canvas)
            self.scan(sec="LI19", loc="244", axis="Y", plot_axes=self.l3YPlot1_axes, plot_canvas=self.l3YPlot1_canvas)
            self.scan(sec="LI19", loc="344", axis="Y", plot_axes=self.l3YPlot1_axes, plot_canvas=self.l3YPlot1_canvas)
        elif area == "IP" and axis == "XY":
            self.singlesscan(axis="X", progressBar=self.progressBar_1)
            self.singlesscan(axis="Y", progressBar=self.progressBar_1)

    def update_quad_k_array(self):
        """Make the quad scan array from k min to k max in n steps."""
        try:
            k_min = float(self.QuadKMinBox.text())
            k_max = float(self.QuadKMaxBox.text())
            n_steps = int(self.QuadStepsBox.text())
        except ValueError:
            return  # a box is empty or not a number yet

        if n_steps < 1:
            return

        self.quad_k_array = np.linspace(k_min, k_max, n_steps)
        print("Quad k array:", self.quad_k_array)

    # This selects thet current active wire.
    def wire_selection_changed(self, selected_wire):
        self.current_wire = selected_wire
        print(f"Selected wire: {selected_wire}")

    def elog_plot(self):
        if self.last_scan is None:
            print("No scan available for Elog.")
            return
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        # outfile = Path(f"wire_scan_{timestamp}.png")
        thescan = self.last_scan
        # logbook = "facet"
        logbook = "mcc"
        Loguser = "Newer Wire Scan Gui"
        LogTitle = f"Wire Scan: {thescan['wire']} ({thescan['axis']})"
        entry_text = (f"Wire: {thescan['wire']}\n"
                      f"Axis: {thescan['axis']}\n"
                      f"Raw RMS: {thescan['raw_rms']:.3f} ± {thescan['raw_rms_error']:.3f} um\n"
                      f"Fitted RMS: {thescan['fitted_rms']:.3f} ± {thescan['fitted_rms_error']:.3f} um\n"
                      f"Skew: {thescan['fitted_skew']:.3f} \n"
                      f"Kurtosis: {thescan['fitted_kurtosis']:.3f}\n"
                      f"Time Stamp: {timestamp}"
                      )
        attachment = str(self.outfile)
        submit_entry(logbook, Loguser, LogTitle, entry_text, attachment)
        print(f"Submitted {logbook} at {attachment} to the Logbook.")

    def calculate_emittance_from_rms(self, section="L2", axis="X"):
        section = section.upper()
        axis = axis.upper()
        section_wires = {
            "L2": ["WIRE:LI11:444", "WIRE:LI11:614", "WIRE:LI11:744", "WIRE:LI12:214"],
            "L3": ["WIRE:LI18:944", "WIRE:LI19:144", "WIRE:LI19:244", "WIRE:LI19:344"],
        }

        axis_index = self.rms_axis_indices[axis]
        wire_indices = [self.rms_wire_indices[wire] for wire in section_wires[section]]
        fitted_rms = self.RMSdata[axis_index, wire_indices, 1]
        fitted_rms_error = self.RMSdata[axis_index, wire_indices, 3]
        # self.fitted_rms = fitted_rms
        # self.fitted_rms_error = fitted_rms_error

        if not np.all(np.isfinite(fitted_rms)):
            raise ValueError(f"Run all four {section} {axis} scans before calculating emittance.")

        # emittance_calculator = EmittClass(fitted_rms.reshape(4, 1))
        emittance_calculator = EmittClass(M=fitted_rms, area=section, axis=axis, )

        beam_twiss = emittance_calculator.EmittIT()
        propagated = emittance_calculator.EmittError(fitted_rms_error)

        emittance = propagated["emittance"]
        emittance_error = propagated["emittance_error"]
        match = propagated["bmag"]
        match_error = propagated["bmag_error"]
        matched_emittance = match * emittance
        matched_emittance_variance = (match ** 2 * emittance_error ** 2 + emittance ** 2 * match_error ** 2
                                      + 2.0 * emittance * match * propagated["covariance_emittance_bmag"])
        matched_emittance_error = np.sqrt(max(matched_emittance_variance, 0.0))

        results = {
            "section": section,
            "axis": axis,
            "fitted_rms": fitted_rms.copy(),
            "fitted_rms_error": fitted_rms_error.copy(),
            "beam_twiss": beam_twiss.copy(),
            "emittance": emittance,
            "emittance_error": emittance_error,
            "match": match,
            "match_error": match_error,
            "matched_emittance": matched_emittance,
            "matched_emittance_error": matched_emittance_error,
        }
        self.emittance_results[(section, axis)] = results
        return results

    def run_emittance_calculation(self, section="L2", axis="X"):
        section = section.upper()
        axis = axis.upper()

        field_suffix = "" if axis == "X" else "_2"
        emittance_field = getattr(self, f"{section.lower()}EmittanceValue{field_suffix}")
        matching_field = getattr(self, f"{section.lower()}MatchingValue{field_suffix}")
        matched_field = getattr(self, f"{section.lower()}MatchedEmittanceValue{field_suffix}")

        try:
            results = self.calculate_emittance_from_rms(section, axis)
        except (KeyError, ValueError) as error:
            message = f"{axis}: {error}"
            emittance_field.setText(message)
            matching_field.clear()
            matched_field.clear()
            print(message)
            return None

        # EmittClass returns geometric emittance in m·rad. Convert it to
        # mm·mrad for display in the existing GUI fields.
        display_scale = 1e6
        emittance_field.setText(
            f"{axis}: {results['emittance'] * display_scale:.4g} ± "
            f"{results['emittance_error'] * display_scale:.3g} mm·mrad"
        )
        matching_field.setText(f"{axis}: {results['match']:.4g} ± " f"{results['match_error']:.3g}")
        matched_field.setText(
            f"{axis}: {results['matched_emittance'] * display_scale:.4g} ± "
            f"{results['matched_emittance_error'] * display_scale:.3g} mm·mrad"
        )
        return results

    def parkFunc(self, axis):
        _, sec, loc = self.WireDropdown.currentText().split(":")
        wireFunk(wirePV=f"WIRE:{sec}:{loc}", axis=axis).ParkerFunk()

    def PARKXY(self):
        if self.axis is None:
            return
        self.parkFunc(axis=self.axis)


class wireFunk:
    def __init__(self, wirePV, area=None, axis="X", shotnum=200, scanvel=100.0, t_out=60.0):
        self.wirePV = wirePV

        # self.pmtPV = pmtPV
        self.dev, self.sec, self.loc = wirePV.split(":")
        self.axis = axis.upper()

        if area is None:
            if self.sec in ("LI11", "LI12"):
                self.area = "L2"
            elif self.sec in ("LI18", "LI19"):
                self.area = "L3"
            else:
                self.area = "IP"
        else:
            self.area = area.upper()

        if self.axis not in ("X", "Y"):
            raise ValueError("axis must be 'X' or 'Y'")
        if not isinstance(shotnum, (int, np.integer)) or shotnum < 3:
            raise ValueError("shotnum must be an integer of at least 3")
        if not np.isfinite(scanvel) or scanvel <= 0:
            raise ValueError("scanvel must be positive")
        if not np.isfinite(t_out) or t_out <= 0:
            raise ValueError("t_out must be positive")

        self.original_scan_limits = None
        self.shotnum = shotnum
        self.scanvel = scanvel
        self.t_out = t_out
        self.pv_x = get_pv(f"{self.wirePV}:USEXWIRE")
        self.pv_y = get_pv(f"{self.wirePV}:USEYWIRE")
        self.OGupper_pv = get_pv(f"{self.wirePV}:{self.axis}WIREOUTER")
        self.OGlower_pv = get_pv(f"{self.wirePV}:{self.axis}WIREINNER")
        self.pv_rate = get_pv("EVNT:SYS1:1:BEAMRATE")
        self.motor_pos = get_pv(f"{self.wirePV}:POSNHSTBR")
        self.lvdt_pos = get_pv(f"{self.wirePV}:LVPOSHSTBR")
        self.pv_startscan = get_pv(f"{self.wirePV}:STARTSCAN")
        self.scanstat_pv = get_pv(f"{self.wirePV}:SCANPROGRESS")
        self.park_pv = get_pv(f"{self.wirePV}:{self.axis}WIREPARKPOSN")
        self.motor_set_pv = get_pv(f"{self.wirePV}:MOTR")
        self.motor_vel_pv = get_pv(f"{self.wirePV}:MOTR.VELO")

        # Selecting PMT PV configuration:
        if self.sec in ("LI11", "LI12"):
            gateW_pvSet = "QADC:LI11:100:TWID"
            if self.loc == "444" or self.loc == "614":
                gate_val = 151
                self.pv_pmt = get_pv("PMT:LI11:614:QDCRAWHSTBR")
                caput(gateW_pvSet, gate_val)
            # elif
            # gate_val = 200
            # self.pv_pmt = get_pv("PMT:LI11:614:QDCRAWHSTBR")
            # caput(gateW_pvSet, gate_val)
            elif self.loc == "744":
                gate_val = 200
                self.pv_pmt = get_pv("PMT:LI11:614:QDCRAWHSTBR")
                caput(gateW_pvSet, gate_val)

            elif self.loc == "214":
                gate_val = 200
                self.pv_pmt = get_pv("PMT:LI11:614:QDCRAWHSTBR")

        elif self.sec in ("LI18", "LI19"):
            self.pv_pmt = get_pv(f"PMT:{self.sec}:{self.loc}:QDCRAWHSTBR")
            gateW_pvSet = "QADC:LI19:100:TWID"

            if self.loc == "944":
                gate_val = 190
                caput(gateW_pvSet, gate_val)

            elif self.loc == "144":
                gate_val = 200
                caput(gateW_pvSet, gate_val)

            elif self.loc == "244":
                gate_val = 100
                caput(gateW_pvSet, gate_val)

            elif self.loc == "344":
                gate_val = 150
                caput(gateW_pvSet, gate_val)

        elif self.sec == "LI20":
            gate_val = 100
            self.pv_pmt = get_pv("PMT:LI20:3350:QDCRAWHSTBR")
            gateW_pvSet = "QADC:LI20:100:TWID"
            caput(gateW_pvSet, gate_val)

    # def scanfunk(self, reverse=None):
    def scanfunk(self, reverse=None, progress_callback=None):
        rate = self.pv_rate.get()
        if rate is None or not np.isfinite(rate) or rate <= 0:
            raise ValueError("Beam rate must be positive")

        if reverse is None or reverse == 0:
            pass
        elif reverse == 1:
            self.reverseit()
        else:
            raise ValueError("Reverse must be 1 to reverse wire scan, or 0 or none for a regular scan")

        self.pv_x.put(int(self.axis == "X"), wait=True)
        self.pv_y.put(int(self.axis == "Y"), wait=True)

        # Stating the Scan
        self.pv_startscan.put(1, wait=True)
        print(f"OK Now Starting {self.axis} wirescan...")
        
        deadline = monotonic() + self.t_out

        while True:
            scan_progress = self.scanstat_pv.get()
        
            if scan_progress is not None:
                scan_progress = max(0, min(100, int(scan_progress)))
        
                if progress_callback is not None:
                    abort_requested = progress_callback(scan_progress)
                    if abort_requested:
                        print("Scan aborted. Collecting available data.")
                        break
        
                if scan_progress >= 100:
                    break
        
            sleep(0.01)

            if monotonic() >= deadline:
                print("Scan timed out. Collecting available data.")
                break

        sleep(0.1)

        if reverse == 1:
            self.restoreit()

        # Keep every PMT sample in the selected history; replace NaN/inf with zero.
        history = self.pv_pmt.get(use_monitor=False)
        history = np.asarray([] if history is None else history, dtype=float).reshape(-1)
        history = np.nan_to_num(history, nan=0.0, posinf=0.0, neginf=0.0)
        data_length = min(self.shotnum, history.size)
        PMTdata = history[-data_length:].copy() if data_length else np.zeros(0)

        # Relative position [µm] = sample time [s] * scan speed [µm/s].
        # scanvel defaults to 100 µm/s; no extra 1000 or 1e6 conversion.
        # Motor/LVDT histories do not limit this time-based position array.
        x = np.arange(data_length, dtype=float) * self.scanvel / rate / np.sqrt(2.0)

        print(f"Completed {self.axis} wirescan.")
        return x, PMTdata #, Motordata #, LVDTdata

    def ParkerFunk(self):
        park_position = self.park_pv.get()
        self.motor_vel_pv.put(1000, wait=True)
        self.motor_set_pv.put(park_position, wait=True)
        print(f"Parked {self.wirePV} {self.axis}")

    def reverseit(self):
        original_outer = self.OGupper_pv.get()
        original_inner = self.OGlower_pv.get()
        # Save the original values only once.
        if self.original_scan_limits is None:
            self.original_scan_limits = (original_outer, original_inner,)
        self.OGupper_pv.put(original_inner, wait=True)
        self.OGlower_pv.put(original_outer, wait=True)

        print(f"Reversed {self.area} {self.axis} scan direction")

    def restoreit(self):
        if self.original_scan_limits is None:
            return
        original_outer, original_inner = self.original_scan_limits

        self.OGupper_pv.put(original_outer, wait=True)
        self.OGlower_pv.put(original_inner, wait=True)

        self.original_scan_limits = None

        print(f"Restored {self.area} {self.axis} scan direction")


class PF_FINDER:
    """Analyze a scan in three steps: detect peaks, calculate RMS, plot.

    Optional calculate_single_scan_error() adds bootstrap uncertainty.
    Public analysis methods come first, shared calculation helpers follow,
    and the plotting entry points and renderer are at the end.
    """

    def __init__(self, x, y, wirePV, wire_diameter_um=45.0):
        self.wirePV = wirePV
        self.x = np.asarray(x, dtype=float)
        self.x = np.nan_to_num(self.x, nan=0.0, posinf=0.0, neginf=0.0)

        self.y = np.asarray(y, dtype=float)
        self.y = np.nan_to_num(self.y, nan=0.0, posinf=0.0, neginf=0.0)

        if self.x.ndim != 1 or self.y.ndim != 1:
            raise ValueError("x and y must be one-dimensional arrays.")
        if len(self.x) != len(self.y):
            raise ValueError("x and y must have the same length.")
        if not np.all(np.isfinite(self.x)) or not np.all(np.isfinite(self.y)):
            raise ValueError("x and y cannot contain NaN or infinite values.")

        self.wire_diameter_um = float(wire_diameter_um)
        if not np.isfinite(self.wire_diameter_um) or self.wire_diameter_um < 0:
            raise ValueError("wire_diameter_um must be finite and nonnegative.")
        self.wire_rms_um = self.wire_diameter_um / 4.0

        self.results = None
        self.peak_windows = None
        self.detection_info = None
        self.beam_result = None
        self.beam_error_result = None
        self.fitter = self  # Compatibility for existing finder.fitter calls.

    # Public analysis

    def find_peaks_and_widths(
            self, prominence=None, height=None, distance=None, rel_height=0.9, noise_factor=5.0,
            min_prominence_fraction=0.05, ):
        """Find physically prominent peaks and their approximate boundaries.

        When ``prominence`` is omitted, a robust threshold is selected from
        the point-to-point noise and the signal's percentile range. Supply a
        numeric ``prominence`` to override the automatic threshold.
        """
        if len(self.x) < 3:
            self.results = []
            self.peak_windows = np.empty((0, 3), dtype=int)
            return self.results, self.peak_windows

        if prominence is None:
            prominence_used, noise_sigma, robust_range = self._automatic_prominence(
                noise_factor=noise_factor, min_prominence_fraction=min_prominence_fraction
            )
            prominence_mode = "automatic"
        else:
            prominence_used = prominence
            noise_sigma = np.nan
            robust_range = np.percentile(self.y, 95) - np.percentile(self.y, 5)
            prominence_mode = "manual"

        self.detection_info = {
            "prominence_mode": prominence_mode,
            "prominence_used": float(prominence_used),
            "estimated_noise_sigma": float(noise_sigma),
            "robust_signal_range": float(robust_range),
            "noise_factor": float(noise_factor),
            "min_prominence_fraction": float(min_prominence_fraction),
        }

        peaks, properties = find_peaks(self.y, prominence=prominence_used, height=height, distance=distance)

        if len(peaks) == 0:
            self.results = []
            self.peak_windows = np.empty((0, 3), dtype=int)
            return self.results, self.peak_windows

        widths, width_heights, left_ips, right_ips = peak_widths(self.y, peaks, rel_height=rel_height)

        sample_indices = np.arange(len(self.x))
        left_x = np.interp(left_ips, sample_indices, self.x)
        right_x = np.interp(right_ips, sample_indices, self.x)

        results = []
        prominences = properties.get("prominences", np.full(len(peaks), np.nan))

        for i, peak_index in enumerate(peaks):
            left_index = max(int(np.floor(left_ips[i])), 0)
            right_index = min(int(np.ceil(right_ips[i])), len(self.x) - 1)

            results.append({
                "peak_number": i + 1,
                "peak_index": int(peak_index),
                "left_index": left_index,
                "right_index": right_index,
                "left_index_fractional": left_ips[i],
                "right_index_fractional": right_ips[i],
                "x_peak": self.x[peak_index],
                "y_peak": self.y[peak_index],
                "x_left": left_x[i],
                "x_right": right_x[i],
                "x_width": right_x[i] - left_x[i],
                "width_samples": widths[i],
                "width_height": width_heights[i],
                "prominence": prominences[i],
            })

        self.results = results
        self.peak_windows = np.array([[peak["left_index"], peak["peak_index"], peak["right_index"]] for peak in self.results], dtype=int)

        return self.results, self.peak_windows

    def calculate_beam_rms(self, raw_y, peak_index=None, smoothed_y=None, fit_asym_gauss=True):
        """Calculate measured and wire-corrected RMS sizes (x in micrometers).

        Existing beam_rms/raw_rms/fitted_rms keys remain measured values.
        corrected_raw_rms/corrected_fitted_rms subtract (wire diameter / 4)
        in quadrature. A zero result with *_wire_resolved=False means the
        measured RMS does not exceed the wire response.

        ``raw_y`` is used for the background estimate and RMS weights.  By
        default, the signal supplied when constructing the class is used as
        the smoothed signal for reporting and plotting. When
        ``fit_asym_gauss`` is true, an asymmetric Gaussian is fitted
        independently to the complete, untrimmed ``x`` and ``raw_y`` arrays.
        Its RMS, R-squared, and residual RMSE are added to the returned result.
        When ``peak_index`` is omitted, the detected peak with the greatest
        prominence is selected automatically. Pass an integer to select a
        different detected peak explicitly.
        """
        if self.results is None:
            raise RuntimeError("Run find_peaks_and_widths() before calculate_beam_rms().")
        if len(self.results) == 0:
            raw_y, smoothed_y = self._prepare_signals(raw_y, smoothed_y)
            # No usable peak: show zero results without inventing extra samples.
            self.beam_result = dict.fromkeys(("baseline", "beam_centroid", "beam_rms", "raw_rms", "fitted_rms", "fitted_centroid", "beam_skew", "beam_kurtosis", "fitted_skew", "fitted_kurtosis", "goodness_of_fit", "corrected_raw_rms", "corrected_fitted_rms"), 0.0)
            self.beam_result.update({"left_x": float(self.x[0]) if self.x.size else 0.0, "right_x": float(self.x[-1]) if self.x.size > 1 else 1.0, "raw_y": raw_y, "smoothed_y": smoothed_y, "fit_result": None, "valid_measurement": False, "beam_data": pd.DataFrame({"x": self.x, "raw_signal": raw_y, "smoothed_signal": smoothed_y})})
            return self.beam_result
        if peak_index is None:
            peak_index = int(np.nanargmax([peak["prominence"] for peak in self.results]))
        if not 0 <= peak_index < len(self.results):
            raise IndexError("peak_index is outside the detected peak range.")

        raw_y, smoothed_y = self._prepare_signals(raw_y, smoothed_y)

        peak = self.results[peak_index]
        left_x = peak["x_left"]
        right_x = peak["x_right"]
        beam_mask = (self.x >= left_x) & (self.x <= right_x)
        background_mask = ~beam_mask

        if not np.any(beam_mask):
            raise ValueError("No samples fall inside the selected beam bounds.")
        if not np.any(background_mask):
            raise ValueError("No background samples remain outside the beam bounds.")

        beam_x = self.x[beam_mask]
        beam_y_raw = raw_y[beam_mask]
        beam_y_smooth = smoothed_y[beam_mask]
        baseline = float(np.median(raw_y[background_mask]))
        beam_weights, beam_centroid, beam_rms = self._weighted_moments(beam_x, beam_y_raw, baseline,
                                                                       "No positive beam signal remains after baseline subtraction.", )
        # SKEW AND KURT SECTION
        centered_x = beam_x - beam_centroid
        total_weight = beam_weights.sum()

        third_moment = np.sum(
            beam_weights * centered_x ** 3
        ) / total_weight

        fourth_moment = np.sum(
            beam_weights * centered_x ** 4
        ) / total_weight

        if beam_rms > 0:
            beam_skew = third_moment / beam_rms ** 3
            beam_kurtosis = fourth_moment / beam_rms ** 4 - 3.0
        else:
            beam_skew = 0.0
            beam_kurtosis = 0.0
        ######

        beam_data = pd.DataFrame(
            {"x": beam_x, "raw_signal": beam_y_raw, "smoothed_signal": beam_y_smooth, "weight": beam_weights})

        fit_result = None
        fitted_rms = 0.0
        fitted_centroid = 0.0
        fitted_skew = 0.0
        fitted_kurtosis = 0.0
        goodness_of_fit = 0.0
        fit_rmse = 0.0

        if fit_asym_gauss and len(self.x) >= 5:
            # Keep the fitted interpretation independent of PF_FINDER's
            # selected beam window: the fitter receives the original full scan.
            fit_result = self.fit_asym_gauss(self.x, raw_y)
            fitted_signal = np.asarray(fit_result["fit_y"], dtype=float)
            fitted_weights, fitted_centroid, fitted_rms = self._weighted_moments(
                self.x,
                fitted_signal,
                fit_result["offset"],
                "No positive fitted beam signal remains after fitted offset subtraction.",
            )

            fitted_centered_x = self.x - fitted_centroid
            fitted_total_weight = fitted_weights.sum()

            fitted_third_moment = np.sum(
                fitted_weights * fitted_centered_x ** 3
            ) / fitted_total_weight

            fitted_fourth_moment = np.sum(
                fitted_weights * fitted_centered_x ** 4
            ) / fitted_total_weight

            if fitted_rms > 0:
                fitted_skew = fitted_third_moment / fitted_rms ** 3
                fitted_kurtosis = (
                        fitted_fourth_moment / fitted_rms ** 4 - 3.0
                )
            # _, fitted_centroid, fitted_rms = self._weighted_moments(
            #     self.x, fitted_signal, fit_result["offset"],
            #     "No positive fitted beam signal remains after fitted offset subtraction.",
            # )

            residuals = raw_y - fitted_signal
            residual_sum_squares = float(np.sum(residuals ** 2))
            total_sum_squares = float(np.sum((raw_y - np.mean(raw_y)) ** 2))
            goodness_of_fit = float(1.0 - residual_sum_squares / total_sum_squares) if total_sum_squares > 0 else 0.0
            fit_rmse = float(np.sqrt(np.mean(residuals ** 2)))
            fit_result["r_squared"] = goodness_of_fit
            fit_result["fit_rmse"] = fit_rmse
            fit_result["fitted_rms"] = fitted_rms
            fit_result["fitted_centroid"] = fitted_centroid
            fit_result["fitted_skew"] = fitted_skew
            fit_result["fitted_kurtosis"] = fitted_kurtosis
        # self.fitted_rms = fitted_rms
        # Keep measured keys unchanged; corrected keys remove the wire response.
        self.beam_result = {
            "peak_index": peak_index,
            "peak": peak,
            "left_x": float(left_x),
            "right_x": float(right_x),
            "baseline": baseline,
            "beam_centroid": beam_centroid,
            "beam_rms": beam_rms,
            "raw_rms": beam_rms,
            "beam_skew": float(beam_skew),  #### Recently Added
            "beam_kurtosis": float(beam_kurtosis),  #### Recently Added
            "fitted_skew": float(fitted_skew),  #### Recently Added
            "fitted_kurtosis": float(fitted_kurtosis),  #### Recently Added
            "fitted_centroid": fitted_centroid,
            "fitted_rms": fitted_rms,
            "goodness_of_fit": goodness_of_fit,
            "r_squared": goodness_of_fit,
            "fit_rmse": fit_rmse,
            "fit_result": fit_result,
            "beam_mask": beam_mask,
            "beam_data": beam_data,
            "raw_y": raw_y,
            "smoothed_y": smoothed_y,
        }
        self.beam_result.update({
            "wire_diameter_um": self.wire_diameter_um,
            "wire_rms_um": self.wire_rms_um,
            "corrected_raw_rms": float(self._correct_wire_rms(beam_rms)),
            "corrected_fitted_rms": float(self._correct_wire_rms(fitted_rms)),
            "raw_wire_resolved": bool(beam_rms > self.wire_rms_um),
            "fitted_wire_resolved": bool(fitted_rms > self.wire_rms_um),
        })
        return self.beam_result

    def calculate_single_scan_error(
            self, raw_y, smoothed_y=None, peak_index=None, n_bootstrap=1000, confidence=0.95, random_state=None
    ):
        """Estimate RMS uncertainty for one scan using residual bootstrap.

        Noise is estimated from ``raw_y - smoothed_y``. Each bootstrap trial
        resamples those centered residuals, adds them to the smoothed scan,
        and recalculates both the windowed raw RMS and the independent
        full-scan asymmetric-Gaussian fitted RMS.

        The returned ``raw_rms_error`` and ``fitted_rms_error`` are the sample
        standard deviations of the successful bootstrap results (one-sigma
        uncertainty). Percentile confidence intervals are also returned.
        """
        if not isinstance(n_bootstrap, (int, np.integer)) or n_bootstrap < 2:
            raise ValueError("n_bootstrap must be an integer of at least 2.")
        if not 0 < confidence < 1:
            raise ValueError("confidence must be between 0 and 1.")

        raw_y, smoothed_y = self._prepare_signals(raw_y, smoothed_y)

        nominal_result = self.calculate_beam_rms(
            raw_y=raw_y, smoothed_y=smoothed_y, peak_index=peak_index, fit_asym_gauss=True
        )
        if not nominal_result.get("valid_measurement", True):
            self.beam_error_result = {"raw_rms_error": 0.0, "fitted_rms_error": 0.0, "corrected_raw_rms_error": 0.0, "corrected_fitted_rms_error": 0.0}
            return self.beam_error_result
        selected_peak_index = nominal_result["peak_index"]
        beam_mask = nominal_result["beam_mask"]
        background_mask = ~beam_mask

        residuals = raw_y - smoothed_y
        centered_residuals = residuals - np.mean(residuals)
        noise_sigma = float(np.std(centered_residuals, ddof=1))
        rng = np.random.default_rng(random_state)
        raw_rms_samples = []
        fitted_rms_samples = []

        for _ in range(n_bootstrap):
            sampled_noise = rng.choice(centered_residuals, size=centered_residuals.size, replace=True)
            synthetic_raw = smoothed_y + sampled_noise

            baseline = float(np.median(synthetic_raw[background_mask]))
            try:
                _, raw_centroid, raw_rms = self._weighted_moments(
                    self.x[beam_mask], synthetic_raw[beam_mask], baseline,
                )
            except ValueError:
                continue

            try:
                fit_result = self.fit_asym_gauss(self.x, synthetic_raw)
            except (RuntimeError, ValueError, FloatingPointError):
                continue

            try:
                _, fitted_centroid, fitted_rms = self._weighted_moments(
                    self.x, np.asarray(fit_result["fit_y"]), fit_result["offset"],
                )
            except ValueError:
                continue
            raw_rms_samples.append(raw_rms)
            fitted_rms_samples.append(fitted_rms)

        successful_trials = len(raw_rms_samples)
        if successful_trials < 2:
            raise RuntimeError("Fewer than two bootstrap trials succeeded; uncertainty " "cannot be estimated.")

        raw_rms_samples = np.asarray(raw_rms_samples)
        fitted_rms_samples = np.asarray(fitted_rms_samples)
        alpha = 1.0 - confidence
        percentiles = [100.0 * alpha / 2.0, 100.0 * (1.0 - alpha / 2.0)]

        error_result = {
            "peak_index": selected_peak_index,
            "n_bootstrap_requested": n_bootstrap,
            "n_bootstrap_successful": successful_trials,
            "confidence": float(confidence),
            "estimated_noise_sigma": noise_sigma,
            "raw_rms": nominal_result["raw_rms"],
            "raw_rms_error": float(np.std(raw_rms_samples, ddof=1)),
            "raw_rms_confidence_interval": tuple(np.percentile(raw_rms_samples, percentiles)),
            "fitted_rms": nominal_result["fitted_rms"],
            "fitted_rms_error": float(np.std(fitted_rms_samples, ddof=1)),
            "fitted_rms_confidence_interval": tuple(np.percentile(fitted_rms_samples, percentiles)),
            "raw_rms_samples": raw_rms_samples,
            "fitted_rms_samples": fitted_rms_samples,
        }
        # Transform every bootstrap draw: the correction is nonlinear.
        error_result["wire_rms_um"] = self.wire_rms_um
        for name, samples in (("raw", raw_rms_samples), ("fitted", fitted_rms_samples)):
            corrected = self._correct_wire_rms(samples)
            prefix = f"corrected_{name}_rms"
            error_result[prefix] = nominal_result[prefix]
            error_result[prefix + "_samples"] = corrected
            error_result[prefix + "_error"] = float(np.std(corrected, ddof=1))
            error_result[prefix + "_confidence_interval"] = tuple(np.percentile(corrected, percentiles))
            error_result[prefix + "_unresolved_trials"] = int(np.count_nonzero(samples <= self.wire_rms_um))
        self.beam_error_result = error_result
        return error_result

    def fit_asym_gauss(self, x, y):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)

        if x.ndim != 1 or y.ndim != 1:
            raise ValueError("x and y must be one-dimensional arrays.")
        if len(x) != len(y):
            raise ValueError("x and y must have the same length.")
        if len(x) < 5:
            raise ValueError("At least five points are recommended for fitting.")
        if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
            raise ValueError("x and y cannot contain NaN or infinite values.")

        peak_i = int(np.argmax(y))
        x_peak = x[peak_i]
        y_peak = y[peak_i]
        offset0 = np.min(y)
        x_range = np.max(x) - np.min(x)
        if x_range <= 0:
            raise ValueError("x must span a nonzero range.")

        p0 = [y_peak - offset0, x_peak, x_range / 6.0, x_range / 6.0, offset0]
        lower_bounds = [0, np.min(x), np.finfo(float).eps, np.finfo(float).eps, -np.inf]
        upper_bounds = [np.inf, np.max(x), x_range, x_range, np.inf]
        popt, pcov = curve_fit(self.asym_gauss, x, y, p0=p0, bounds=(lower_bounds, upper_bounds), maxfev=10000)

        amp, x0, sig_left, sig_right, offset = popt
        effective_sigma = np.sqrt((sig_left ** 2 + sig_right ** 2) / 2)
        fwhm = np.sqrt(2 * np.log(2)) * (sig_left + sig_right)
        y_fit = self.asym_gauss(x, *popt)
        curve_start_x, curve_end_x = self.curve_bounds_from_params(x0, sig_left, sig_right)
        return {
            "peak_index_local": peak_i,
            "peak_x_raw": x_peak,
            "peak_y_raw": y_peak,
            "fit_amp": amp,
            "fit_peak_x": x0,
            "sigma_left": sig_left,
            "sigma_right": sig_right,
            "effective_sigma": effective_sigma,
            "full_width_half_max": fwhm,
            "offset": offset,
            "curve_threshold_fraction": 0.05,
            "curve_start_x": curve_start_x,
            "curve_end_x": curve_end_x,
            "fit_x": x,
            "fit_y": y_fit,
            "fit_params": popt,
            "fit_covariance": pcov,
        }

    # Shared calculations
    def _automatic_prominence(self, noise_factor=5.0, min_prominence_fraction=0.05):
        """Estimate a prominence that rejects baseline-scale fluctuations."""
        if noise_factor <= 0:
            raise ValueError("noise_factor must be positive.")
        if not 0 < min_prominence_fraction < 1:
            raise ValueError("min_prominence_fraction must be between 0 and 1.")

        # First differences remove a slowly varying baseline. Dividing by
        # sqrt(2) converts their robust spread to per-sample noise sigma.
        differences = np.diff(self.y)
        difference_median = np.median(differences)
        difference_mad = np.median(np.abs(differences - difference_median))
        noise_sigma = 1.4826 * difference_mad / np.sqrt(2.0)

        # Percentiles keep a single bad sample from defining the full range.
        robust_range = np.percentile(self.y, 95) - np.percentile(self.y, 5)
        noise_prominence = noise_factor * noise_sigma
        relative_prominence = min_prominence_fraction * robust_range
        automatic_prominence = max(noise_prominence, relative_prominence)

        return float(automatic_prominence), float(noise_sigma), float(robust_range)

    def _correct_wire_rms(self, measured_rms):
        """Remove (diameter / 4) in quadrature; x must be in micrometers.

        Negative corrected variance is clipped to zero (wire-unresolved),
        not interpreted as a resolved zero-sized beam. NaN becomes zero.
        The wire diameter is treated as exact in the bootstrap uncertainty.
        """
        # WIRE SIZE SUBTRACTED IN QUADRATURE HERE: sqrt(measured² - (diameter/4)²).
        # This uses the existing round-wire response model and 45 µm default diameter.
        measured_rms = np.nan_to_num(measured_rms, nan=0.0, posinf=0.0, neginf=0.0)
        return np.sqrt(np.maximum(np.asarray(measured_rms) ** 2 - self.wire_rms_um ** 2, 0.0))

    def _prepare_signals(self, raw_y, smoothed_y):
        """Normalize and validate signals for nominal and bootstrap analysis."""
        raw_y = np.nan_to_num(np.asarray(raw_y, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
        if smoothed_y is None:
            smoothed_y = self.y
        else:
            smoothed_y = np.nan_to_num(np.asarray(smoothed_y, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)

        if raw_y.shape != self.x.shape or smoothed_y.shape != self.x.shape:
            raise ValueError("raw_y and smoothed_y must have the same shape as x.")
        if not np.all(np.isfinite(raw_y)) or not np.all(np.isfinite(smoothed_y)):
            raise ValueError("raw_y and smoothed_y cannot contain NaN or infinite values.")

        return raw_y, smoothed_y

    @staticmethod
    def _weighted_moments(x, signal, baseline,
                          empty_message="No positive signal remains after baseline subtraction."):
        """Return clipped weights, centroid and RMS using sample-weighted sums.

        Both the selected raw window and full fitted scan use this formula;
        bootstrap trials use exactly the same calculation.
        """
        weights = np.clip(signal - baseline, 0, None)
        total_weight = weights.sum()
        if total_weight <= 0:
            raise ValueError(empty_message)
        centroid = float(np.sum(x * weights) / total_weight)
        rms = float(np.sqrt(np.sum(weights * (x - centroid) ** 2) / total_weight))
        return weights, centroid, rms

    @staticmethod
    def asym_gauss(x, amp, x0, sig_left, sig_right, offset):
        sigma = np.where(x < x0, sig_left, sig_right)
        return offset + amp * np.exp(-0.5 * ((x - x0) / sigma) ** 2)

    @staticmethod
    def curve_bounds_from_params(x0, sig_left, sig_right, threshold_fraction=0.05):
        if not 0 < threshold_fraction < 1:
            raise ValueError("threshold_fraction must be between 0 and 1.")
        if sig_left <= 0 or sig_right <= 0:
            raise ValueError("sig_left and sig_right must be positive.")

        multiplier = np.sqrt(-2.0 * np.log(threshold_fraction))
        return (float(x0 - multiplier * sig_left), float(x0 + multiplier * sig_right))

    def plot_beam_rms(self, beam_result=None, error_result=None, figsize=(10, 5), show=True, ax=None):
        """Plot RMS and bootstrap uncertainty through the shared renderer.

        The signature and (figure, axes) return value are preserved. Omitted
        results use the most recently calculated beam and bootstrap results.
        """
        # return self.plot(mode="beam", beam_result=beam_result, error_result=error_result, figsize=figsize, show=show, ax=ax)
        return self.plot(beam_result=beam_result, error_result=error_result, figsize=figsize, show=show, ax=ax, )

    def plot(self, beam_result=None, error_result=None,
             figsize=(10, 5), show=True, ax=None):
        """Renderer for the beam view.

        All artist creation, styling, labels and presentation live here.
        Existing axes are cleared and reused, including embedded Qt axes.
        """
        if beam_result is None:
            beam_result = self.beam_result
        if beam_result is None:
            raise RuntimeError("Run calculate_beam_rms() before plotting, or pass a beam_result.")

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.figure
            ax.clear()

        beam_data = beam_result["beam_data"]
        beam_rms = beam_result["beam_rms"]
        beam_centroid = beam_result["beam_centroid"]
        baseline = beam_result["baseline"]
        fitted_rms = beam_result.get("fitted_rms", 0.0)
        beam_skew = beam_result.get("beam_skew", 0.0)
        beam_kurtosis = beam_result.get("beam_kurtosis", 0.0)
        fitted_skew = beam_result.get("fitted_skew", 0.0)
        fitted_kurtosis = beam_result.get("fitted_kurtosis", 0.0)
        goodness_of_fit = beam_result.get("goodness_of_fit", np.nan)
        plot_x = self.x if "raw_y" in beam_result else beam_data["x"]

        if error_result is None:
            error_result = self.beam_error_result

        raw_rms_error = np.nan
        fitted_rms_error = np.nan

        if error_result is not None:
            raw_rms_error = error_result.get("raw_rms_error", np.nan)
            fitted_rms_error = error_result.get("fitted_rms_error", np.nan)

        ax.scatter(plot_x, beam_result.get("raw_y", beam_data["raw_signal"]), s=20, label="Raw beam data")
        # Moving AVG
        # ax.plot(plot_x, beam_result.get("smoothed_y", beam_data["smoothed_signal"]), color="orange", linewidth=2,label="Moving-average", )

        fit_result = beam_result.get("fit_result")
        if fit_result is not None:
            ax.plot(fit_result["fit_x"], fit_result["fit_y"], color="deepskyblue", linewidth=2.5,
                    label="Asym Gaussian fit", )

        ax.axvspan(beam_result["left_x"], beam_result["right_x"], color="mediumpurple", alpha=0.15,
                   label="RMS Region", )

        # Not Complex - Text Box Constructor
        ######################################
        result_lines = [
            f"Raw RMS = {beam_rms:.3f} ± {raw_rms_error:.3f} µm",
            f"Fitted RMS = {fitted_rms:.3f} ± {fitted_rms_error:.3f} µm",
            f"Raw skew/kurtosis = {beam_skew:.3f} / {beam_kurtosis:.3f}",
            f"Fitted skew/kurtosis = {fitted_skew:.3f} / {fitted_kurtosis:.3f}",
        ]
        # ax.text(0.02, 0.98, "\n".join(result_lines), transform=ax.transAxes, va="top", fontsize=9,bbox={"facecolor": "white", "alpha": 0.8}, )
        ######################################

        ax.axhline(baseline, color="gray", linestyle="--", label="Baseline")
        ax.axvline(beam_centroid, color="limegreen", linestyle="--", label=f"Centroid = {beam_centroid:.3g}")

        # Title Construction:
        ######################################
        title = (f" {self.wirePV} \n"
                 f" RMS {fitted_rms:.3f} ±{fitted_rms_error:.3f} | Skew: {beam_skew:.3f} | kurtosis: {beam_kurtosis:.3f} ")
        title += (f"\nWire-corrected raw/fit RMS: {beam_result['corrected_raw_rms']:.3f} / {beam_result['corrected_fitted_rms']:.3f} µm")
        if not beam_result.get("valid_measurement", True):
            title += "\nNo usable peak (zero placeholders)"
        self.title = title
        # Configuring the axis and Title.
        ax.set_title(title)
        ax.set_xlabel("Relative wire travel (µm)")
        ax.set_ylabel("Cnts")
        ######################################

        peaks = [beam_result["peak"]] if "peak" in beam_result else []

        for peak in peaks:
            ax.plot(peak["x_peak"], peak["y_peak"], "x", markersize=10)
            ax.axvline(peak["x_left"], linestyle=":", linewidth=2)
            ax.axvline(peak["x_right"], linestyle=":", linewidth=2)
            ax.hlines(y=peak["width_height"], xmin=peak["x_left"],
                      xmax=peak["x_right"], linewidth=2)
            # ax.annotate(f"Peak {peak['peak_number']}", xy=(peak["x_peak"], peak["y_peak"]), xytext=(0, 10), textcoords="offset points", ha="center")
        left_x = beam_result["left_x"]
        right_x = beam_result["right_x"]
        padding = 0.5 * (right_x - left_x)

        ax.set_xlim(left_x - padding, right_x + padding, )
        # ax.legend(loc="best")
        fig.tight_layout()
        fig.canvas.draw_idle()

        if show:
            plt.show()
        return fig, ax


class EmittClass:
    def __init__(self, M, area, axis):
        self.M = np.asarray(M, dtype=float).reshape(-1)
        if self.M.size != 4:
            raise ValueError("M must contain four beam-size measurements")
        self.area = area.upper()
        self.axis = axis.upper()

        self.wire_names = ["WS11444", "WS11614", "WS11744", "WS12214",
                           "WS18944", "WS19144", "WS19244", "WS19344", "IPWS1", ]
        if self.area == "L2":
            k = 0
        elif self.area == "L3":
            k = 4
        else:
            raise ValueError("area must be 'L2' or 'L3'")

        if self.axis == "X":
            j = 0
        elif self.axis == "Y":
            j = 1
        else:
            raise ValueError("axis must be 'X' or 'Y'")

        # Columns: Z [m], Beta_X [m], Beta_Y [m], Alpha_X, Alpha_Y, Eta_X, Eta_Y , Phase Adv [deg].
        self.DTwiss = np.array([
            [1065.0098, 7.892, 33.41, 0.7829, -3.798, -1.380e-07, 0.0, 0],  # L2
            [1080.2526, 34.49, 16.88, 1.178, -1.024, -7.141e-08, 0.0, 90],  # L2
            [1102.0430, 46.06, 28.87, -1.642, 0.5384, 1.377e-07, 0.0, 180],  # L2
            [1132.4750, 66.87, 26.98, 2.057, -0.589, 2.538e-07, 0.0, 270],  # L2
            [1826.9139, 11.44, 46.65, -0.4504, 2.351, -2.032e-06, 0.0, 0],  # L3
            [1840.9015, 44.61, 8.26, -1.921, 0.3942, 6.927e-06, 0.0, 90],  # L3
            [1853.2459, 12.46, 37.00, 0.6501, -2.023, 7.249e-06, 0.0, 180],  # L3
            [1865.5903, 42.72, 11.83, -1.917, 0.4439, 1.618e-05, 0.0, 270],  # L3
            [1993.9100, 0.5032, 0.5032, -0.08, -0.08, -3.285e-07, -4.118e-15, 0],  # IPWS1
        ], dtype=float)
        # Columns: Gamma_X[m], Gamma_Y [m]
        self.gammas = np.column_stack((
            (1 + self.DTwiss[:, 3] ** 2) / self.DTwiss[:, 1],  # Gamma X
            (1 + self.DTwiss[:, 4] ** 2) / self.DTwiss[:, 2],  # Gamma Y
        ))

        self.Twiss = np.array(
            [
                [self.DTwiss[0 + k, 1 + j], self.DTwiss[1 + k, 1 + j], self.DTwiss[2 + k, 1 + j],
                 self.DTwiss[3 + k, 1 + j]],  # beta
                [self.DTwiss[0 + k, 3 + j], self.DTwiss[1 + k, 3 + j], self.DTwiss[2 + k, 3 + j],
                 self.DTwiss[3 + k, 3 + j]],  # alpha
                [self.DTwiss[0 + k, 7], self.DTwiss[1 + k, 7], self.DTwiss[2 + k, 7], self.DTwiss[3 + k, 7]],
                # mu (Phase advance) in degrees
                [self.gammas[0 + k, 0 + j], self.gammas[1 + k, 0 + j], self.gammas[2 + k, 0 + j],
                 self.gammas[3 + k, 0 + j]],  # gamma = (1 + alpha**2)/beta
            ]
        )

        self.SIG2 = (self.M * 1e-6) ** 2
        self.m11_ = np.empty(4)
        self.m12_ = np.empty(4)
        self.MIO = np.empty((4, 3))

        Beta_0 = self.Twiss[0, 0]
        Alpha_0 = self.Twiss[1, 0]
        Mu_0 = self.Twiss[2, 0]

        for k in range(4):
            Beta_k = self.Twiss[0, k]
            Mu_k = self.Twiss[2, k]
            Delta_mu = np.deg2rad(Mu_k - Mu_0)

            self.m11_[k] = np.sqrt(Beta_k / Beta_0) * (np.cos(Delta_mu) + Alpha_0 * np.sin(Delta_mu))
            self.m12_[k] = np.sqrt(Beta_0 * Beta_k) * np.sin(Delta_mu)

            self.MIO[k, :] = [self.m11_[k] ** 2, 2 * self.m11_[k] * self.m12_[k], self.m12_[k] ** 2]

        print("MIO:")
        print(self.MIO)

    def EmittIT(self):
        self.A, self.B, self.C = np.linalg.lstsq(self.MIO, self.SIG2, rcond=None)[0]

        determinant = self.A * self.C - self.B ** 2

        if determinant <= 0:
            raise ValueError("A*C - B**2 is not positive. The fitted beam matrix is unphysical.")

        self.Emitt = np.sqrt(determinant)
        self.Beta_beam = self.A / self.Emitt
        self.Alpha_beam = -self.B / self.Emitt
        self.Gamma_beam = self.C / self.Emitt

        Beta_0 = self.Twiss[0, 0]
        Alpha_0 = self.Twiss[1, 0]
        Gamma_0 = self.Twiss[3, 0]

        self.Bmag = 0.5 * (self.Beta_beam * Gamma_0 - 2 * self.Alpha_beam * Alpha_0 + self.Gamma_beam * Beta_0)

        self.BeamTwiss = np.array([self.Emitt, self.Beta_beam, self.Alpha_beam, self.Gamma_beam, self.Bmag])

        print("Emittance:", self.Emitt)
        print("Bmag:", self.Bmag)

        return self.BeamTwiss

    def EmittError(self, M_error):
        """Return propagated 1-sigma errors for emittance and Bmag."""
        M_error = np.asarray(M_error, dtype=float)
        if M_error.shape != (4,):
            raise ValueError("M_error must contain one uncertainty for each of the " "four measurements.")
        if not hasattr(self, "Emitt"):
            raise RuntimeError("Run EmittIT() before EmittError().")

        sigma = self.M * 1e-6
        sigma_error = M_error * 1e-6
        sig2_error = 2.0 * np.abs(sigma) * sigma_error
        covariance_sig2 = np.diag(sig2_error ** 2)

        fit_operator = np.linalg.pinv(self.MIO)
        covariance_abc = fit_operator @ covariance_sig2 @ fit_operator.T

        A, B, C = self.A, self.B, self.C
        emittance = self.Emitt
        grad_emittance = np.array([C / (2.0 * emittance), -B / emittance, A / (2.0 * emittance)])

        beta_0 = self.Twiss[0, 0]
        alpha_0 = self.Twiss[1, 0]
        gamma_0 = self.Twiss[3, 0]
        numerator_gradient = np.array([gamma_0, 2.0 * alpha_0, beta_0])
        grad_bmag = numerator_gradient / (2.0 * emittance) - (self.Bmag / emittance) * grad_emittance

        emittance_variance = grad_emittance @ covariance_abc @ grad_emittance
        bmag_variance = grad_bmag @ covariance_abc @ grad_bmag
        emittance_bmag_covariance = grad_emittance @ covariance_abc @ grad_bmag

        return {
            "emittance": emittance,
            "emittance_error": np.sqrt(max(emittance_variance, 0.0)),
            "bmag": self.Bmag,
            "bmag_error": np.sqrt(max(bmag_variance, 0.0)),
            "covariance_emittance_bmag": emittance_bmag_covariance,
            "covariance_ABC": covariance_abc,
        }


def main():
    app = QtWidgets.QApplication(sys.argv)
    window = MyGui()
    window.show()
    return app.exec_()


if __name__ == "__main__":
    sys.exit(main())

