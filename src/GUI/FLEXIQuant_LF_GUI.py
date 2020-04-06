#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FLEXIQuant LF GUI version:
This script creates a graphical user interface for FLEXIQuant LF using PyQT5
and connects it to the FLEXIQuant_LF_method_GUI.py script
"""

# Author: Konstantin Kahnert
# Date: 2019_10_15
# Python version: 3.7.3

# # import required libraries

import sys
# used to build the GUI
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import *
from PyQt5.QtCore import pyqtSignal
# used to check is file path exists and to adapt it's format
from os.path import exists, normpath
# used to execute the FLEXIQuant LF method script as a new process
from subprocess import Popen
# FLEXIQuant LF method script
from src.GUI.FLEXIQuant_LF_method_GUI import *
# used for regular expressions
from re import findall
# used to save plots as pdf document
from matplotlib import use
use("pdf")


# # define classes

# define worker thread class to execute the FLEXIQuant LF method script in
class FQLFWorkerThread(QtCore.QThread):

    # initiate signals to communicate between main and worker thread
    sig_progress = pyqtSignal(float)
    sig_run_button = pyqtSignal(bool)
    sig_cancel_button = pyqtSignal(bool)
    sig_error_ref_samples = pyqtSignal(str)
    sig_error_file_open = pyqtSignal(str)
    sig_error_group_column = pyqtSignal(str)

    # define init function
    def __init__(self, ref_samples, path_input, filename_input, path_output,
                 num_ransac_init, mod_cutoff, make_plots, parent=None):
        QtCore.QThread.__init__(self, parent)
        self.ref_samples = ref_samples
        self.path_input = path_input
        self.filename_input = filename_input
        self.path_output = path_output
        self.num_ransac_init = num_ransac_init
        self.mod_cutoff = mod_cutoff
        self.make_plots = make_plots

    def send_progress(self, progress):
        """Sends progress as signal from worker thread to main thread"""
        self.sig_progress.emit(progress)

    def send_error_ref_samples(self, name):
        """Sends reference sample error signal from worker thread to main thread"""
        self.sig_error_ref_samples.emit(name)

    def send_error_file_open(self, path):
        """Sends file open error signal from worker thread to main thread"""
        self.sig_error_file_open.emit(path)

    def send_error_group_column(self):
        """ Sends group column signal from worker thread to main thread"""
        self.sig_error_group_column.emit("")

    def enable_run_button(self):
        """ Sends signal from worker thread to main thread to make the run button clickable again"""
        self.sig_run_button.emit(True)

    def disable_cancel_button(self):
        """ Sends signal from worker thread to main thread to make cancel button unclickable"""
        self.sig_cancel_button.emit(False)

    def stop(self):
        """Terminates the running process and resets progressbar as well as run and cancel buttons"""
        self.send_progress(0)
        self.enable_run_button()
        self.disable_cancel_button()
        self.terminate()

    def run(self):
        """Executes the FLEXIQuant LF method script"""
        run_FQLF_main(self, ui)


# define main window class
class Ui_MainWindow(object):

    # define setup function
    def setupUi(self, MainWindow):
        """Defines how the main window looks and creates it"""
        # set object name
        MainWindow.setObjectName("MainWindow")

        # set fixed window size
        MainWindow.setFixedSize(700, 410)

        # initiate and name main window widget
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")

        # set FQLF logo as window icon
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(r"D:\Masterarbeit\FLEXIQuant_LF_logo.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.setWindowIcon(icon)

        # initiate and name GridLayout widget
        self.gridLayout = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName("gridLayout")

        # initiate label for input file selection field
        self.label_input = QtWidgets.QLabel(self.centralwidget)
        # set font
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.label_input.setFont(font)
        # define text type (rich text)
        self.label_input.setTextFormat(QtCore.Qt.RichText)
        # name object
        self.label_input.setObjectName("label_input")
        # set grid position
        self.gridLayout.addWidget(self.label_input, 0, 0, 1, 1)

        # initiate help button
        self.pushButton_help = QtWidgets.QPushButton(self.centralwidget)
        # set object name
        self.pushButton_help.setObjectName("pushButton")
        # make button clickable
        self.pushButton_help.setEnabled(True)
        # connect to help_button_clicked function (sets what happens when clicked)
        self.pushButton_help.clicked.connect(self.help_button_clicked)
        # set grid position
        self.gridLayout.addWidget(self.pushButton_help, 26, 5, 1, 1)

        # initiate label for reference sample identifier input field
        self.label_ref_samples = QtWidgets.QLabel(self.centralwidget)
        # set font
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.label_ref_samples.setFont(font)
        # define text type (rich text)
        self.label_ref_samples.setTextFormat(QtCore.Qt.RichText)
        # set object name
        self.label_ref_samples.setObjectName("label_ref_samples")
        # set grid position
        self.gridLayout.addWidget(self.label_ref_samples, 5, 0, 1, 1)

        # initiate empty field to enter input file path
        self.lineEdit_output = QtWidgets.QLineEdit(self.centralwidget)
        # set object name
        self.lineEdit_output.setObjectName("lineEdit_output")
        # set grid position
        self.gridLayout.addWidget(self.lineEdit_output, 4, 0, 1, 5)

        # initiate label for output folder path input field
        self.label_output = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.label_output.setFont(font)
        self.label_output.setTextFormat(QtCore.Qt.RichText)
        # set object name
        self.label_output.setObjectName("label_output")
        # set grid position
        self.gridLayout.addWidget(self.label_output, 3, 0, 1, 1)

        # initiate a horizontal spacer item
        spacerItem_horizontal = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding,
                                                      QtWidgets.QSizePolicy.Minimum)
        # set grid position
        self.gridLayout.addItem(spacerItem_horizontal, 12, 3, 1, 1)

        # initiate a spinbox to select number of ransac initiations
        self.spinBox_ransac_init = QtWidgets.QSpinBox(self.centralwidget)
        # set font size
        self.spinBox_ransac_init.setStyleSheet("font: 8pt")
        # set size of widget to a fixed size
        self.spinBox_ransac_init.setMinimumSize(QtCore.QSize(55, 22))
        self.spinBox_ransac_init.setMaximumSize(QtCore.QSize(55, 22))
        # set minimum and maximum value that can be selected
        self.spinBox_ransac_init.setMinimum(5)
        self.spinBox_ransac_init.setMaximum(100)
        # set step size for clicking up or down
        self.spinBox_ransac_init.setSingleStep(5)
        # set default value to 30 initiations
        self.spinBox_ransac_init.setProperty("value", 30)
        # set object name
        self.spinBox_ransac_init.setObjectName("spinBox_ransac_init")
        # set grid position
        self.gridLayout.addWidget(self.spinBox_ransac_init, 12, 1, 1, 1)

        # initiate label for options section
        self.label_options = QtWidgets.QLabel(self.centralwidget)
        # set font
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.label_options.setFont(font)
        # set text type
        self.label_options.setTextFormat(QtCore.Qt.RichText)
        # set object name
        self.label_options.setObjectName("label_options")
        # set grid position
        self.gridLayout.addWidget(self.label_options, 8, 0, 1, 1)

        # initiate empty field to enter reference sample identifier
        self.lineEdit_ref_samples = QtWidgets.QLineEdit(self.centralwidget)
        # set object name
        self.lineEdit_ref_samples.setObjectName("lineEdit_ref_samples")
        # set grid position
        self.gridLayout.addWidget(self.lineEdit_ref_samples, 6, 0, 1, 5)

        # initiate label for github repo link
        self.label_github_link = QtWidgets.QLabel(self.centralwidget)
        # set font
        font = QtGui.QFont()
        font.setPointSize(7)
        font.setBold(False)
        font.setItalic(False)
        font.setUnderline(True)
        font.setWeight(50)
        self.label_github_link.setFont(font)
        # set text type
        self.label_github_link.setTextFormat(QtCore.Qt.RichText)
        # # set maximum size, scaling and word wrap, not needed with fixed size window
        # self.label_github_link.setMaximumSize(QtCore.QSize(16777215, 20))
        # self.label_github_link.setScaledContents(False)
        # self.label_github_link.setWordWrap(False)
        # set object name
        self.label_github_link.setObjectName("label_github_link")
        # add hyperlink to github repo
        self.label_github_link.setOpenExternalLinks(True)
        github_link = "<a href=\"https://github.com/SteenOmicsLab/FLEXIQuantLF\">GitHub Repository</a>"
        self.label_github_link.setText(github_link)
        # set grid position
        self.gridLayout.addWidget(self.label_github_link, 28, 5, 1, 1)

        # initiate progress bar
        self.progressBar = QtWidgets.QProgressBar(self.centralwidget)
        # set default value to 0
        self.progressBar.setProperty("value", 0)
        # set object name
        self.progressBar.setObjectName("progressBar")
        # set grid position
        self.gridLayout.addWidget(self.progressBar, 27, 0, 1, 4)

        # initiate cancel button
        self.pushButton_cancel = QtWidgets.QPushButton(self.centralwidget)
        # set object name
        self.pushButton_cancel.setObjectName("pushButton_cancel")
        # make cancel button unclickable
        self.pushButton_cancel.setEnabled(False)
        # set grid position
        self.gridLayout.addWidget(self.pushButton_cancel, 27, 5, 1, 1)

        # add spacer item between input field and options
        spacerItem_middle = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum,
                                                  QtWidgets.QSizePolicy.Expanding)
        # set grid position
        self.gridLayout.addItem(spacerItem_middle, 7, 0, 1, 1)

        # initiate label for number of ransac initiations spinbox
        self.label_ransac_init = QtWidgets.QLabel(self.centralwidget)
        # set object name
        self.label_ransac_init.setObjectName("label_ransac_init")
        # set grid position
        self.gridLayout.addWidget(self.label_ransac_init, 12, 0, 1, 1, QtCore.Qt.AlignLeft)

        # initiate select output button
        self.pushButton_select_output = QtWidgets.QPushButton(self.centralwidget)
        # set object name
        self.pushButton_select_output.setObjectName("pushButton_select_output")
        # connect to select_output_button_clicked function (sets what happens when clicked)
        self.pushButton_select_output.clicked.connect(self.select_output_button_clicked)
        # set grid position
        self.gridLayout.addWidget(self.pushButton_select_output, 4, 5, 1, 1)

        # initiate run button
        self.pushButton_run = QtWidgets.QPushButton(self.centralwidget)
        # set object name
        self.pushButton_run.setObjectName("pushButton_run")
        # make button clickable
        self.pushButton_run.setEnabled(False)
        # connect to run_button_clicked function (sets what happens when clicked)
        self.pushButton_run.clicked.connect(self.run_button_clicked)
        # set grid position
        self.gridLayout.addWidget(self.pushButton_run, 27, 4, 1, 1)

        # initiate empty field to enter input file path
        self.lineEdit_input = QtWidgets.QLineEdit(self.centralwidget)
        # set object name
        self.lineEdit_input.setObjectName("lineEdit_input")
        # set grid position
        self.gridLayout.addWidget(self.lineEdit_input, 2, 0, 1, 5)

        # initiate select input button
        self.pushButton_select_input = QtWidgets.QPushButton(self.centralwidget)
        # set object name
        self.pushButton_select_input.setObjectName("pushButton_select_input")
        # connect to select_input_button_clicked function (sets what happens when clicked)
        self.pushButton_select_input.clicked.connect(self.select_input_button_clicked)
        # set grid position
        self.gridLayout.addWidget(self.pushButton_select_input, 2, 5, 1, 1)

        # initiate label for modification cutoff spinbox
        self.label_mod_cutoff = QtWidgets.QLabel(self.centralwidget)
        # set object name
        self.label_mod_cutoff.setObjectName("label_mod_cutoff")
        # set grid position
        self.gridLayout.addWidget(self.label_mod_cutoff, 13, 0, 1, 1)

        # initiate spinbox to select modification cutoff
        self.doubleSpinBox_mod_cutoff = QtWidgets.QDoubleSpinBox(self.centralwidget)
        # set font size
        self.doubleSpinBox_mod_cutoff.setStyleSheet("font: 8pt")
        # set size of widget to a fixed size
        self.doubleSpinBox_mod_cutoff.setMinimumSize(QtCore.QSize(55, 22))
        self.doubleSpinBox_mod_cutoff.setMaximumSize(QtCore.QSize(55, 22))
        # set number of decimals displayed
        self.doubleSpinBox_mod_cutoff.setDecimals(2)
        # set maximum value that can be selected
        self.doubleSpinBox_mod_cutoff.setMaximum(1.0)
        # set step size for clicking up and down
        self.doubleSpinBox_mod_cutoff.setSingleStep(0.05)
        # set default value to 0.5
        self.doubleSpinBox_mod_cutoff.setProperty("value", 0.5)
        # set object name
        self.doubleSpinBox_mod_cutoff.setObjectName("doubleSpinBox_mod_cutoff")
        # set grid position
        self.gridLayout.addWidget(self.doubleSpinBox_mod_cutoff, 13, 1, 1, 1)

        # set style sheet for check box
        StyleSheet_checkBox_plots = '''
        QCheckBox::indicator {
            width:  16px;
            height: 16px;
        }
        '''
        # initiate check box to select if plots should be created or not
        self.checkBox_plots = QtWidgets.QCheckBox(self.centralwidget)
        # use style sheet created above
        self.checkBox_plots.setStyleSheet(StyleSheet_checkBox_plots)
        # set fixed size of widget
        self.checkBox_plots.setMinimumSize(QtCore.QSize(16, 16))
        self.checkBox_plots.setMaximumSize(QtCore.QSize(16, 16))
        # set orientation of text and checkbox
        self.checkBox_plots.setLayoutDirection(QtCore.Qt.RightToLeft)
        # set text of check box to nothing
        self.checkBox_plots.setText("")
        # set default state to unchecked
        self.checkBox_plots.setTristate(False)
        # set object name
        self.checkBox_plots.setObjectName("checkBox_plots")
        # set grid position
        self.gridLayout.addWidget(self.checkBox_plots, 14, 1, 1, 1, QtCore.Qt.AlignRight)

        # initiate label for checkbox
        self.label_create_plots = QtWidgets.QLabel(self.centralwidget)
        # set object name
        self.label_create_plots.setObjectName("label_create_plots")
        # set grid position
        self.gridLayout.addWidget(self.label_create_plots, 14, 0, 1, 1)

        # set central widget of mainwindow
        MainWindow.setCentralWidget(self.centralwidget)

        # call retranslateUi function
        self.retranslateUi(MainWindow)

        # connect slots by name
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        """Tries to translate and set the text that will be displayed in the GUI to the .
        If translation is not successful, english source text below will be displayed"""
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "FLEXIQuant LF"))
        self.label_input.setText(_translate("MainWindow", "Input file"))
        self.pushButton_help.setText(_translate("MainWindow", "Help"))
        self.label_ref_samples.setText(_translate("MainWindow", "Reference sample identifier"))
        self.label_output.setText(_translate("MainWindow", "Output folder"))
        self.label_options.setText(_translate("MainWindow", "Options"))
        self.pushButton_cancel.setText(_translate("MainWindow", "Cancel"))
        self.label_ransac_init.setText(_translate("MainWindow", "RANSAC initiations:"))
        self.pushButton_select_output.setText(_translate("MainWindow", "Select"))
        self.pushButton_run.setText(_translate("MainWindow", "Run"))
        self.pushButton_select_input.setText(_translate("MainWindow", "Select"))
        self.label_mod_cutoff.setText(_translate("MainWindow", "Modification cutoff:"))
        self.label_create_plots.setText(_translate("MainWindow", "Create plots:"))


    def select_input_button_clicked(self):
        """ Sets what happens when the select input button is clicked
        Opens a file browser where the user can select the input file (only displays .csv files)
        displays the file path of the selected file in input path field and the folder path in output folder field
        checks if input file and output folder filed are filled and paths are valid
        and if true, makes run button clickable"""
        # initiate QFileDialog widget
        filedialog = QFileDialog()

        # open file browser and save path of selected file
        self.fpath_input = filedialog.getOpenFileName(filedialog, filter="CSV files (*.csv)")[0]

        # display path of select file in lineEdit_input
        self.lineEdit_input.setText(self.fpath_input)

        # get file name
        self.fname_input = findall("[^\\\\,/]*$", self.fpath_input)[0]

        # get input folder
        self.folder_output = self.fpath_input[0:-len(self.fname_input)-1]

        # display input folder in lineEdit_output (as default)
        self.lineEdit_output.setText(self.folder_output)

        # check if input and output fields are filled
        if exists(self.folder_output) is True and exists(self.fpath_input) is True:
            self.pushButton_run.setEnabled(True)


    def select_output_button_clicked(self):
        """ Sets what happens when the select output button is clicked
        Opens a file browser where the user can select the output folder
        displays the file path of the selected file in input path field and the folder path in output folder field
        checks if input file and output folder filed are filled and paths are valid
        and if true, makes run button clickable"""
        # initiate QFileDialog widget
        filedialog = QFileDialog()

        # open file browser and save path of selected file
        self.folder_output = filedialog.getExistingDirectory(filedialog, "Select Directory")

        # display path of select folder in lineEdit_output
        self.lineEdit_output.setText(self.folder_output)

        # check if input and output fields are filled
        if exists(self.folder_output) is True and exists(self.fpath_input) is True:
            self.pushButton_run.setEnabled(True)


    def error_ref_samples(self, name):
        """ Displays error message and resets progressbar and run and cancel buttons"""
        # display error message
        self.display_error_msg("Given reference sample identifier \"" + name + "\" not found in \"Group\" column")
        # reset progressbar and run and cancel buttons
        self.reset_after_run_error()

    def error_file_open(self, path):
        """ Displays error message and resets progressbar and run and cancel buttons"""
        # display error message
        self.display_error_msg("Permission denied!\n" + "Please close " + path)
        # reset progressbar and run and cancel buttons
        self.reset_after_run_error()

    def error_group_column(self, key):
        """ Displays error message and resets progressbar and run and cancel buttons"""
        # display error message
        self.display_error_msg("Incorrect input format!\n" + "No column named \"Group\" found in input file.")
        # reset progressbar and run and cancel buttons
        self.reset_after_run_error()

    def reset_after_run_error(self):
        """ Resets progressbar and run and cancel buttons"""
        # enable start button
        self.update_run_button(True)

        # disable cancel button
        self.pushButton_cancel.setEnabled(False)

        # set progress bar to 0
        self.progressBar.setValue(0)

    def update_progress(self, progress):
        """ Updates the progressbar to the given integer (progress)
        and displays message when program has finished"""
        self.progressBar.setValue(progress)

        if progress == 100:
            self.display_done_msg()

    def update_run_button(self, on_off):
        """Makes run button clickable or unclickable
        Arguments:
        True: clickable
        False: unclickable"""
        self.pushButton_run.setEnabled(on_off)

    def update_cancel_button(self, on_off):
        """Makes cancel button clickable or unclickable
        Arguments:
        True: clickable
        False: unclickable"""
        self.pushButton_cancel.setEnabled(on_off)

    def display_done_msg(self):
        """ Opens dialog window with task completion message and asks user to open output folder """
        # enable start button
        self.update_run_button(True)

        # disable cancel button
        self.pushButton_cancel.setEnabled(False)

        # initiate dialog box
        self.msg = QMessageBox()
        # set message box type
        self.msg.setIcon(QMessageBox.Information)
        # set window title
        self.msg.setWindowTitle("Task completed")
        # set text
        self.msg.setText("Done!")
        # set informative text
        self.msg.setInformativeText("Open output folder?")
        # add yes and no buttons
        self.msg.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        # set yes button as default selection
        self.msg.setDefaultButton(QMessageBox.Yes)

        # display dialog
        reply = self.msg.exec_()

        # if Yes is clicked, open output folder
        if reply == QMessageBox.Yes:
            # open output folder
            output_dir = normpath(self.path_output)
            Popen(r'explorer "' + output_dir + '"')
        else:
            pass

    def display_error_msg(self, msg_txt):
        """ Opens dialog window and displays error message with msg_text
        msg_text: str, error message to be displayed"""
        # initiate message box
        self.msg = QMessageBox()
        # set icon
        self.msg.setIcon(QMessageBox.Critical)
        # set text
        self.msg.setText("Error")
        # set informative text
        self.msg.setInformativeText(msg_txt)
        # set window title
        self.msg.setWindowTitle("Error")
        # display error message window
        self.msg.exec_()

    def help_button_clicked(self):
        """Opens documentation"""
        Popen([r'FLEXIQuantLF_GUI_Documentation.pdf'], shell=True)

    def run_button_clicked(self):
        """ Sets what happens when the run button is clicked
        collects all paths and parameters set by the user and displays error message if paths don't exist
        starts the FLEXIQuant_LF_method script in a new thread with the given parameters
        """
        # get reference samples
        ref_samples = self.lineEdit_ref_samples.text()

        # get path of input file
        path_input = self.lineEdit_input.text()

        # get file name
        filename_input = self.fname_input

        # get path of output folder
        self.path_output = self.lineEdit_output.text()

        # get number of ransac initiations
        num_ransac_init = self.spinBox_ransac_init.value()

        # get modification cutoff
        mod_cutoff = self.doubleSpinBox_mod_cutoff.value()

        # get checkBox_plots value
        if self.checkBox_plots.isChecked() == True:
            make_plots = True
        else:
            make_plots = False

        # check if input and output fields are filled
        if exists(path_input) is True and exists(self.path_output) is True:

            # disable start button
            self.update_run_button(False)

            # enable cancel button
            self.pushButton_cancel.setEnabled(True)

            # initiate worker thread with parameters set by user
            self.worker_thread = FQLFWorkerThread(ref_samples, path_input, filename_input, self.path_output,
                                                  num_ransac_init, mod_cutoff, make_plots)
            # start worker thread
            self.worker_thread.start()

            # connect cancel button to terminate thread function
            self.pushButton_cancel.clicked.connect(self.worker_thread.stop)

            # connect to all signals from worker thread
            self.worker_thread.sig_cancel_button.connect(self.update_cancel_button)
            self.worker_thread.sig_run_button.connect(self.update_run_button)
            self.worker_thread.sig_error_ref_samples.connect(self.error_ref_samples)
            self.worker_thread.sig_error_file_open.connect(self.error_file_open)
            self.worker_thread.sig_error_group_column.connect(self.error_group_column)
            self.worker_thread.sig_progress.connect(self.update_progress)

        # if output path doesn't exist, display corresponding error message
        elif exists(path_input) is True and exists(self.path_output) is False:
            self.display_error_msg("Specified output folder does not exist!")

        # if input path doesn't exist, display corresponding error message
        elif exists(path_input) is False and exists(self.path_output) is True:
            self.display_error_msg("Specified input file does not exist!")

        # if input and output folder don't exist, display corresponding error message
        else:
            self.display_error_msg("Specified input file and output folder do not exist!")


# main script
if __name__ == "__main__":
    # initiate QApplication
    app = QtWidgets.QApplication(sys.argv)
    # initiate main window
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()

    # set up main window
    ui.setupUi(MainWindow)

    # display main window
    MainWindow.show()

    # display exit code
    sys.exit(app.exec_())

