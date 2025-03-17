# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'rd1c.ui'
#
# Created by: PyQt5 UI code generator 5.13.0
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets
###***###
from PyQt5.QtWidgets import QWidget, QPushButton, QLineEdit, QInputDialog, QApplication, QLabel, QFileDialog, QMessageBox
from PyQt5.QtCore import QObject, pyqtSignal
from PyQt5.QtGui import QImage, QPixmap
from QtImageViewer_slt1 import QtImageViewer

import sys
import os
import shutil
import re
import time
import glob
#import logging

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors

from sklearn.cluster import KMeans
from sklearn import metrics
from scipy.spatial.distance import cdist
from sklearn.metrics import silhouette_score

__author__ = "Stephen Trisno and Xinghao Zhang"
__version__ = '0.2.2'

#note, need to install/set up environent for all of the above, plus:
#  1)leidenalg


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1269, 868)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.albl_cell = QtWidgets.QLabel(self.centralwidget)
        self.albl_cell.setGeometry(QtCore.QRect(10, 10, 191, 21))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.albl_cell.setFont(font)
        self.albl_cell.setObjectName("albl_cell")
        self.albl_gene = QtWidgets.QLabel(self.centralwidget)
        self.albl_gene.setGeometry(QtCore.QRect(10, 32, 191, 21))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.albl_gene.setFont(font)
        self.albl_gene.setObjectName("albl_gene")
        self.albl_obsgrp = QtWidgets.QLabel(self.centralwidget)
        self.albl_obsgrp.setGeometry(QtCore.QRect(10, 54, 201, 21))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.albl_obsgrp.setFont(font)
        self.albl_obsgrp.setObjectName("albl_obsgrp")
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(10, 132, 101, 21))
        font = QtGui.QFont()
        font.setPointSize(16)
        font.setUnderline(True)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.lbl_filtcell = QtWidgets.QLabel(self.centralwidget)
        self.lbl_filtcell.setGeometry(QtCore.QRect(10, 159, 261, 21))
        font = QtGui.QFont()
        font.setPointSize(11)
        self.lbl_filtcell.setFont(font)
        self.lbl_filtcell.setObjectName("lbl_filtcell")
        self.lbl_filtgene = QtWidgets.QLabel(self.centralwidget)
        self.lbl_filtgene.setGeometry(QtCore.QRect(10, 189, 279, 21))
        font = QtGui.QFont()
        font.setPointSize(11)
        self.lbl_filtgene.setFont(font)
        self.lbl_filtgene.setObjectName("lbl_filtgene")
        self.txt_filtcell = QtWidgets.QTextEdit(self.centralwidget)
        self.txt_filtcell.setEnabled(False)
        self.txt_filtcell.setGeometry(QtCore.QRect(272, 159, 81, 25))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(9)
        self.txt_filtcell.setFont(font)
        self.txt_filtcell.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_filtcell.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_filtcell.setTabChangesFocus(True)
        self.txt_filtcell.setAcceptRichText(False)
        self.txt_filtcell.setObjectName("txt_filtcell")
        self.txt_filtgene = QtWidgets.QTextEdit(self.centralwidget)
        self.txt_filtgene.setEnabled(False)
        self.txt_filtgene.setGeometry(QtCore.QRect(288, 189, 81, 25))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(9)
        self.txt_filtgene.setFont(font)
        self.txt_filtgene.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_filtgene.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_filtgene.setTabChangesFocus(True)
        self.txt_filtgene.setAcceptRichText(False)
        self.txt_filtgene.setObjectName("txt_filtgene")
        self.bttn_recalc = QtWidgets.QPushButton(self.centralwidget)
        self.bttn_recalc.setEnabled(False)
        self.bttn_recalc.setGeometry(QtCore.QRect(10, 85, 111, 33))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(14)
        self.bttn_recalc.setFont(font)
        self.bttn_recalc.setObjectName("bttn_recalc")
        self.lbl_genepercell = QtWidgets.QLabel(self.centralwidget)
        self.lbl_genepercell.setGeometry(QtCore.QRect(10, 234, 111, 21))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.lbl_genepercell.setFont(font)
        self.lbl_genepercell.setObjectName("lbl_genepercell")
        self.txt_genepercell = QtWidgets.QTextEdit(self.centralwidget)
        self.txt_genepercell.setEnabled(False)
        self.txt_genepercell.setGeometry(QtCore.QRect(124, 234, 81, 25))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(8)
        self.txt_genepercell.setFont(font)
        self.txt_genepercell.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_genepercell.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_genepercell.setTabChangesFocus(True)
        self.txt_genepercell.setAcceptRichText(False)
        self.txt_genepercell.setObjectName("txt_genepercell")
        self.lbl_totalgenect = QtWidgets.QLabel(self.centralwidget)
        self.lbl_totalgenect.setGeometry(QtCore.QRect(10, 264, 130, 21))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.lbl_totalgenect.setFont(font)
        self.lbl_totalgenect.setObjectName("lbl_totalgenect")
        self.txt_totalgenect = QtWidgets.QTextEdit(self.centralwidget)
        self.txt_totalgenect.setEnabled(False)
        self.txt_totalgenect.setGeometry(QtCore.QRect(141, 264, 81, 25))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(8)
        self.txt_totalgenect.setFont(font)
        self.txt_totalgenect.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_totalgenect.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_totalgenect.setTabChangesFocus(True)
        self.txt_totalgenect.setAcceptRichText(False)
        self.txt_totalgenect.setPlaceholderText("")
        self.txt_totalgenect.setObjectName("txt_totalgenect")
        self.txt_mitoperc = QtWidgets.QTextEdit(self.centralwidget)
        self.txt_mitoperc.setEnabled(False)
        self.txt_mitoperc.setGeometry(QtCore.QRect(118, 294, 61, 25))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(8)
        self.txt_mitoperc.setFont(font)
        self.txt_mitoperc.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_mitoperc.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_mitoperc.setTabChangesFocus(True)
        self.txt_mitoperc.setAcceptRichText(False)
        self.txt_mitoperc.setObjectName("txt_mitoperc")
        self.lbl_ERCC = QtWidgets.QLabel(self.centralwidget)
        self.lbl_ERCC.setGeometry(QtCore.QRect(10, 324, 81, 21))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.lbl_ERCC.setFont(font)
        self.lbl_ERCC.setObjectName("lbl_ERCC")
        self.lbl_mitoperc = QtWidgets.QLabel(self.centralwidget)
        self.lbl_mitoperc.setGeometry(QtCore.QRect(10, 294, 104, 21))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.lbl_mitoperc.setFont(font)
        self.lbl_mitoperc.setObjectName("lbl_mitoperc")
        self.txt_ERCC = QtWidgets.QTextEdit(self.centralwidget)
        self.txt_ERCC.setEnabled(False)
        self.txt_ERCC.setGeometry(QtCore.QRect(93, 324, 61, 25))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(8)
        self.txt_ERCC.setFont(font)
        self.txt_ERCC.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_ERCC.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_ERCC.setTabChangesFocus(True)
        self.txt_ERCC.setAcceptRichText(False)
        self.txt_ERCC.setObjectName("txt_ERCC")
        self.checkBox_stdgenename = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_stdgenename.setEnabled(False)
        self.checkBox_stdgenename.setGeometry(QtCore.QRect(250, 235, 181, 17))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.checkBox_stdgenename.setFont(font)
        self.checkBox_stdgenename.setChecked(True)
        self.checkBox_stdgenename.setObjectName("checkBox_stdgenename")
        self.checkBox_excludegene = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_excludegene.setEnabled(False)
        self.checkBox_excludegene.setGeometry(QtCore.QRect(250, 256, 191, 61))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.checkBox_excludegene.setFont(font)
        self.checkBox_excludegene.setChecked(False)
        self.checkBox_excludegene.setObjectName("checkBox_excludegene")
        self.checkBox_cellcyclescore = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_cellcyclescore.setEnabled(False)
        self.checkBox_cellcyclescore.setGeometry(QtCore.QRect(250, 320, 181, 17))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.checkBox_cellcyclescore.setFont(font)
        self.checkBox_cellcyclescore.setChecked(True)
        self.checkBox_cellcyclescore.setObjectName("checkBox_cellcyclescore")
        self.lbl_numhivargenes2keep = QtWidgets.QLabel(self.centralwidget)
        self.lbl_numhivargenes2keep.setGeometry(QtCore.QRect(10, 368, 101, 34))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.lbl_numhivargenes2keep.setFont(font)
        self.lbl_numhivargenes2keep.setObjectName("lbl_numhivargenes2keep")
        self.txt_numhivargenes2keep = QtWidgets.QTextEdit(self.centralwidget)
        self.txt_numhivargenes2keep.setEnabled(False)
        self.txt_numhivargenes2keep.setGeometry(QtCore.QRect(115, 375, 71, 25))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(8)
        self.txt_numhivargenes2keep.setFont(font)
        self.txt_numhivargenes2keep.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_numhivargenes2keep.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_numhivargenes2keep.setTabChangesFocus(True)
        self.txt_numhivargenes2keep.setAcceptRichText(False)
        self.txt_numhivargenes2keep.setObjectName("txt_numhivargenes2keep")
        self.lbl_mmmean_vargenes = QtWidgets.QLabel(self.centralwidget)
        self.lbl_mmmean_vargenes.setGeometry(QtCore.QRect(10, 410, 124, 34))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.lbl_mmmean_vargenes.setFont(font)
        self.lbl_mmmean_vargenes.setObjectName("lbl_mmmean_vargenes")
        self.txt_minmean_vargenes = QtWidgets.QTextEdit(self.centralwidget)
        self.txt_minmean_vargenes.setEnabled(False)
        self.txt_minmean_vargenes.setGeometry(QtCore.QRect(144, 415, 61, 25))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(8)
        self.txt_minmean_vargenes.setFont(font)
        self.txt_minmean_vargenes.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_minmean_vargenes.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_minmean_vargenes.setTabChangesFocus(True)
        self.txt_minmean_vargenes.setAcceptRichText(False)
        self.txt_minmean_vargenes.setObjectName("txt_minmean_vargenes")
        self.txt_maxmean_vargenes = QtWidgets.QTextEdit(self.centralwidget)
        self.txt_maxmean_vargenes.setEnabled(False)
        self.txt_maxmean_vargenes.setGeometry(QtCore.QRect(220, 415, 61, 25))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(8)
        self.txt_maxmean_vargenes.setFont(font)
        self.txt_maxmean_vargenes.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_maxmean_vargenes.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_maxmean_vargenes.setTabChangesFocus(True)
        self.txt_maxmean_vargenes.setAcceptRichText(False)
        self.txt_maxmean_vargenes.setObjectName("txt_maxmean_vargenes")
        self.lbl_mmmean_vargenes_slash = QtWidgets.QLabel(self.centralwidget)
        self.lbl_mmmean_vargenes_slash.setGeometry(QtCore.QRect(210, 415, 21, 21))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.lbl_mmmean_vargenes_slash.setFont(font)
        self.lbl_mmmean_vargenes_slash.setObjectName("lbl_mmmean_vargenes_slash")
        self.lbl_mdispersion_vargenes = QtWidgets.QLabel(self.centralwidget)
        self.lbl_mdispersion_vargenes.setGeometry(QtCore.QRect(10, 445, 124, 34))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.lbl_mdispersion_vargenes.setFont(font)
        self.lbl_mdispersion_vargenes.setObjectName("lbl_mdispersion_vargenes")
        self.txt_mdispersion_vargenes = QtWidgets.QTextEdit(self.centralwidget)
        self.txt_mdispersion_vargenes.setEnabled(False)
        self.txt_mdispersion_vargenes.setGeometry(QtCore.QRect(144, 450, 61, 25))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(8)
        self.txt_mdispersion_vargenes.setFont(font)
        self.txt_mdispersion_vargenes.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_mdispersion_vargenes.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_mdispersion_vargenes.setTabChangesFocus(True)
        self.txt_mdispersion_vargenes.setAcceptRichText(False)
        self.txt_mdispersion_vargenes.setObjectName("txt_mdispersion_vargenes")
        self.checkBox_reg_phase = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_reg_phase.setEnabled(False)
        self.checkBox_reg_phase.setGeometry(QtCore.QRect(320, 400, 151, 17))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.checkBox_reg_phase.setFont(font)
        self.checkBox_reg_phase.setCheckable(True)
        self.checkBox_reg_phase.setChecked(False)
        self.checkBox_reg_phase.setObjectName("checkBox_reg_phase")
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(320, 375, 101, 21))
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setUnderline(True)
        self.label_2.setFont(font)
        self.label_2.setObjectName("label_2")
        self.checkBox_reg_mito_nct = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_reg_mito_nct.setEnabled(False)
        self.checkBox_reg_mito_nct.setGeometry(QtCore.QRect(320, 420, 141, 17))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.checkBox_reg_mito_nct.setFont(font)
        self.checkBox_reg_mito_nct.setChecked(False)
        self.checkBox_reg_mito_nct.setObjectName("checkBox_reg_mito_nct")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        self.bttn_next = QtWidgets.QPushButton(self.centralwidget)
        self.bttn_next.setEnabled(False)
        self.bttn_next.setGeometry(QtCore.QRect(130, 85, 88, 33))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        self.bttn_next.setFont(font)
        self.bttn_next.setObjectName("bttn_next")
        self.lbl_NC_listgenes = QtWidgets.QLabel(self.centralwidget)
        self.lbl_NC_listgenes.setGeometry(QtCore.QRect(10, 525, 121, 31))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.lbl_NC_listgenes.setFont(font)
        self.lbl_NC_listgenes.setObjectName("lbl_NC_listgenes")
        self.label_3 = QtWidgets.QLabel(self.centralwidget)
        self.label_3.setGeometry(QtCore.QRect(10, 497, 171, 21))
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setUnderline(True)
        self.label_3.setFont(font)
        self.label_3.setObjectName("label_3")
        self.bttn_NC_listgenes = QtWidgets.QPushButton(self.centralwidget)
        self.bttn_NC_listgenes.setEnabled(False)
        self.bttn_NC_listgenes.setGeometry(QtCore.QRect(130, 525, 71, 33))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        self.bttn_NC_listgenes.setFont(font)
        self.bttn_NC_listgenes.setObjectName("bttn_NC_listgenes")
        self.lbl_NC_neighborhoodnum = QtWidgets.QLabel(self.centralwidget)
        self.lbl_NC_neighborhoodnum.setGeometry(QtCore.QRect(10, 559, 101, 31))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.lbl_NC_neighborhoodnum.setFont(font)
        self.lbl_NC_neighborhoodnum.setObjectName("lbl_NC_neighborhoodnum")
        self.txt_NC_neighborhoodnum = QtWidgets.QTextEdit(self.centralwidget)
        self.txt_NC_neighborhoodnum.setEnabled(False)
        self.txt_NC_neighborhoodnum.setGeometry(QtCore.QRect(108, 563, 33, 25))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(8)
        self.txt_NC_neighborhoodnum.setFont(font)
        self.txt_NC_neighborhoodnum.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_NC_neighborhoodnum.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_NC_neighborhoodnum.setTabChangesFocus(True)
        self.txt_NC_neighborhoodnum.setAcceptRichText(False)
        self.txt_NC_neighborhoodnum.setObjectName("txt_NC_neighborhoodnum")
        self.txt_NC_PCA = QtWidgets.QTextEdit(self.centralwidget)
        self.txt_NC_PCA.setEnabled(False)
        self.txt_NC_PCA.setGeometry(QtCore.QRect(196, 563, 33, 25))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(8)
        self.txt_NC_PCA.setFont(font)
        self.txt_NC_PCA.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_NC_PCA.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_NC_PCA.setTabChangesFocus(True)
        self.txt_NC_PCA.setAcceptRichText(False)
        self.txt_NC_PCA.setObjectName("txt_NC_PCA")
        self.lbl_NC_PCA = QtWidgets.QLabel(self.centralwidget)
        self.lbl_NC_PCA.setGeometry(QtCore.QRect(157, 559, 41, 31))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.lbl_NC_PCA.setFont(font)
        self.lbl_NC_PCA.setObjectName("lbl_NC_PCA")
        self.lbl_NC_resclust = QtWidgets.QLabel(self.centralwidget)
        self.lbl_NC_resclust.setGeometry(QtCore.QRect(10, 592, 151, 36))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.lbl_NC_resclust.setFont(font)
        self.lbl_NC_resclust.setObjectName("lbl_NC_resclust")
        self.txt_NC_resclust = QtWidgets.QTextEdit(self.centralwidget)
        self.txt_NC_resclust.setEnabled(False)
        self.txt_NC_resclust.setGeometry(QtCore.QRect(163, 598, 41, 25))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(8)
        self.txt_NC_resclust.setFont(font)
        self.txt_NC_resclust.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_NC_resclust.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_NC_resclust.setTabChangesFocus(True)
        self.txt_NC_resclust.setAcceptRichText(False)
        self.txt_NC_resclust.setObjectName("txt_NC_resclust")
        self.lbl_NC_dist = QtWidgets.QLabel(self.centralwidget)
        self.lbl_NC_dist.setGeometry(QtCore.QRect(157, 631, 51, 31))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.lbl_NC_dist.setFont(font)
        self.lbl_NC_dist.setObjectName("lbl_NC_dist")
        self.txt_NC_dist = QtWidgets.QTextEdit(self.centralwidget)
        self.txt_NC_dist.setEnabled(False)
        self.txt_NC_dist.setGeometry(QtCore.QRect(210, 635, 41, 25))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(8)
        self.txt_NC_dist.setFont(font)
        self.txt_NC_dist.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_NC_dist.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_NC_dist.setTabChangesFocus(True)
        self.txt_NC_dist.setAcceptRichText(False)
        self.txt_NC_dist.setObjectName("txt_NC_dist")
        self.txt_NC_dimemb = QtWidgets.QTextEdit(self.centralwidget)
        self.txt_NC_dimemb.setEnabled(False)
        self.txt_NC_dimemb.setGeometry(QtCore.QRect(109, 635, 33, 25))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(8)
        self.txt_NC_dimemb.setFont(font)
        self.txt_NC_dimemb.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_NC_dimemb.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_NC_dimemb.setTabChangesFocus(True)
        self.txt_NC_dimemb.setAcceptRichText(False)
        self.txt_NC_dimemb.setObjectName("txt_NC_dimemb")
        self.lbl_NC_dimemb = QtWidgets.QLabel(self.centralwidget)
        self.lbl_NC_dimemb.setGeometry(QtCore.QRect(10, 630, 101, 31))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.lbl_NC_dimemb.setFont(font)
        self.lbl_NC_dimemb.setObjectName("lbl_NC_dimemb")
        self.label_4 = QtWidgets.QLabel(self.centralwidget)
        self.label_4.setGeometry(QtCore.QRect(320, 477, 151, 21))
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setUnderline(True)
        self.label_4.setFont(font)
        self.label_4.setObjectName("label_4")
        self.bttn_cellcyclescore = QtWidgets.QPushButton(self.centralwidget)
        self.bttn_cellcyclescore.setEnabled(False)
        self.bttn_cellcyclescore.setGeometry(QtCore.QRect(440, 314, 71, 33))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        self.bttn_cellcyclescore.setFont(font)
        self.bttn_cellcyclescore.setObjectName("bttn_cellcyclescore")
        self.bttn_excludegene = QtWidgets.QPushButton(self.centralwidget)
        self.bttn_excludegene.setEnabled(False)
        self.bttn_excludegene.setGeometry(QtCore.QRect(440, 270, 71, 33))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        self.bttn_excludegene.setFont(font)
        self.bttn_excludegene.setObjectName("bttn_excludegene")
        self.comboBox_rankgene = QtWidgets.QComboBox(self.centralwidget)
        self.comboBox_rankgene.setEnabled(False)
        self.comboBox_rankgene.setGeometry(QtCore.QRect(320, 510, 151, 25))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.comboBox_rankgene.setFont(font)
        self.comboBox_rankgene.setCurrentText("")
        self.comboBox_rankgene.setObjectName("comboBox_rankgene")
        self.lbl_rankgene = QtWidgets.QLabel(self.centralwidget)
        self.lbl_rankgene.setGeometry(QtCore.QRect(317, 544, 131, 36))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.lbl_rankgene.setFont(font)
        self.lbl_rankgene.setObjectName("lbl_rankgene")
        self.txt_rankgenenum = QtWidgets.QTextEdit(self.centralwidget)
        self.txt_rankgenenum.setEnabled(False)
        self.txt_rankgenenum.setGeometry(QtCore.QRect(440, 550, 33, 25))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(8)
        self.txt_rankgenenum.setFont(font)
        self.txt_rankgenenum.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_rankgenenum.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_rankgenenum.setTabChangesFocus(True)
        self.txt_rankgenenum.setAcceptRichText(False)
        self.txt_rankgenenum.setPlaceholderText("")
        self.txt_rankgenenum.setObjectName("txt_rankgenenum")
        self.textEdit_genevis = QtWidgets.QTextEdit(self.centralwidget)
        self.textEdit_genevis.setEnabled(False)
        self.textEdit_genevis.setGeometry(QtCore.QRect(150, 710, 161, 111))
        self.textEdit_genevis.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.textEdit_genevis.setTabChangesFocus(True)
        self.textEdit_genevis.setAcceptRichText(False)
        self.textEdit_genevis.setObjectName("textEdit_genevis")
        self.lbl_genevis = QtWidgets.QLabel(self.centralwidget)
        self.lbl_genevis.setGeometry(QtCore.QRect(150, 680, 161, 31))
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setUnderline(True)
        self.lbl_genevis.setFont(font)
        self.lbl_genevis.setObjectName("lbl_genevis")
        self.label_5 = QtWidgets.QLabel(self.centralwidget)
        self.label_5.setGeometry(QtCore.QRect(330, 606, 151, 21))
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setUnderline(True)
        self.label_5.setFont(font)
        self.label_5.setObjectName("label_5")
        self.lbl_numgeneplot = QtWidgets.QLabel(self.centralwidget)
        self.lbl_numgeneplot.setGeometry(QtCore.QRect(327, 630, 131, 51))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.lbl_numgeneplot.setFont(font)
        self.lbl_numgeneplot.setObjectName("lbl_numgeneplot")
        self.txt_numgeneplot = QtWidgets.QTextEdit(self.centralwidget)
        self.txt_numgeneplot.setEnabled(False)
        self.txt_numgeneplot.setGeometry(QtCore.QRect(450, 645, 51, 25))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(8)
        self.txt_numgeneplot.setFont(font)
        self.txt_numgeneplot.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_numgeneplot.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txt_numgeneplot.setTabChangesFocus(True)
        self.txt_numgeneplot.setAcceptRichText(False)
        self.txt_numgeneplot.setObjectName("txt_numgeneplot")
        self.bttn_plot = QtWidgets.QPushButton(self.centralwidget)
        self.bttn_plot.setEnabled(False)
        self.bttn_plot.setGeometry(QtCore.QRect(330, 772, 111, 51))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.bttn_plot.setFont(font)
        self.bttn_plot.setObjectName("bttn_plot")
        self.lbl_renameclust = QtWidgets.QLabel(self.centralwidget)
        self.lbl_renameclust.setGeometry(QtCore.QRect(10, 680, 131, 31))
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setUnderline(True)
        self.lbl_renameclust.setFont(font)
        self.lbl_renameclust.setObjectName("lbl_renameclust")
        self.textEdit_renameclust = QtWidgets.QTextEdit(self.centralwidget)
        self.textEdit_renameclust.setEnabled(False)
        self.textEdit_renameclust.setGeometry(QtCore.QRect(10, 710, 131, 111))
        self.textEdit_renameclust.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.textEdit_renameclust.setTabChangesFocus(True)
        self.textEdit_renameclust.setAcceptRichText(False)
        self.textEdit_renameclust.setObjectName("textEdit_renameclust")
        self.checkBox_plot_dot = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_plot_dot.setEnabled(False)
        self.checkBox_plot_dot.setGeometry(QtCore.QRect(330, 687, 71, 17))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.checkBox_plot_dot.setFont(font)
        self.checkBox_plot_dot.setChecked(True)
        self.checkBox_plot_dot.setObjectName("checkBox_plot_dot")
        self.checkBox_plot_stkviolinH = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_plot_stkviolinH.setEnabled(False)
        self.checkBox_plot_stkviolinH.setGeometry(QtCore.QRect(330, 708, 151, 17))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.checkBox_plot_stkviolinH.setFont(font)
        self.checkBox_plot_stkviolinH.setChecked(True)
        self.checkBox_plot_stkviolinH.setObjectName("checkBox_plot_stkviolinH")
        self.checkBox_plot_stkviolinV = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_plot_stkviolinV.setEnabled(False)
        self.checkBox_plot_stkviolinV.setGeometry(QtCore.QRect(330, 729, 151, 17))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.checkBox_plot_stkviolinV.setFont(font)
        self.checkBox_plot_stkviolinV.setChecked(True)
        self.checkBox_plot_stkviolinV.setObjectName("checkBox_plot_stkviolinV")
        self.checkBox_plot_save = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_plot_save.setEnabled(False)
        self.checkBox_plot_save.setGeometry(QtCore.QRect(450, 782, 91, 31))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.checkBox_plot_save.setFont(font)
        self.checkBox_plot_save.setChecked(True)
        self.checkBox_plot_save.setObjectName("checkBox_plot_save")
        self.checkBox_plot_heat = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_plot_heat.setEnabled(False)
        self.checkBox_plot_heat.setGeometry(QtCore.QRect(420, 687, 71, 17))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.checkBox_plot_heat.setFont(font)
        self.checkBox_plot_heat.setChecked(True)
        self.checkBox_plot_heat.setObjectName("checkBox_plot_heat")
        self.checkBox_plot_matrix = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_plot_matrix.setEnabled(False)
        self.checkBox_plot_matrix.setGeometry(QtCore.QRect(330, 750, 81, 17))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.checkBox_plot_matrix.setFont(font)
        self.checkBox_plot_matrix.setChecked(True)
        self.checkBox_plot_matrix.setObjectName("checkBox_plot_matrix")
        self.checkBox_plot_tracks = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_plot_tracks.setEnabled(False)
        self.checkBox_plot_tracks.setGeometry(QtCore.QRect(420, 750, 91, 17))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.checkBox_plot_tracks.setFont(font)
        self.checkBox_plot_tracks.setChecked(True)
        self.checkBox_plot_tracks.setObjectName("checkBox_plot_tracks")
        self.plainTextEdit = QtWidgets.QPlainTextEdit(self.centralwidget)
        self.plainTextEdit.setGeometry(QtCore.QRect(530, 640, 721, 181))
        palette = QtGui.QPalette()
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Light, brush)
        brush = QtGui.QBrush(QtGui.QColor(72, 226, 214))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Text, brush)
        brush = QtGui.QBrush(QtGui.QColor(27, 27, 27))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Light, brush)
        brush = QtGui.QBrush(QtGui.QColor(72, 226, 214))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Text, brush)
        brush = QtGui.QBrush(QtGui.QColor(27, 27, 27))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Light, brush)
        brush = QtGui.QBrush(QtGui.QColor(120, 120, 120))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Text, brush)
        brush = QtGui.QBrush(QtGui.QColor(240, 240, 240))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Base, brush)
        self.plainTextEdit.setPalette(palette)
        font = QtGui.QFont()
        font.setFamily("Courier New")
        font.setPointSize(8)
        self.plainTextEdit.setFont(font)
        self.plainTextEdit.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.plainTextEdit.setReadOnly(True)
        self.plainTextEdit.setObjectName("plainTextEdit")
        self.scrollArea = QtWidgets.QScrollArea(self.centralwidget)
        self.scrollArea.setGeometry(QtCore.QRect(530, 40, 721, 536))
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setObjectName("scrollArea")
        ###***###
        # self.scrollAreaWidgetContents = QtWidgets.QWidget()
        # self.scrollAreaWidgetContents.setGeometry(QtCore.QRect(0, 0, 719, 534))
        # self.scrollAreaWidgetContents.setObjectName("scrollAreaWidgetContents")
        # self.scrollArea.setWidget(self.scrollAreaWidgetContents)
        self.viewer = QtImageViewer()
        self.viewer.setGeometry(QtCore.QRect(0, 0, 719, 534))
        self.viewer.canPan = True
        self.viewer.canZoom = True
        self.viewer.aspectRatioMode = QtCore.Qt.KeepAspectRatio
        # self.viewer.setImage(QPixmap.fromImage(QImage("new_file1.png")))
        self.scrollArea.setWidget(self.viewer)  # scrollAreaWidgetContents)
        MainWindow.setCentralWidget(self.centralwidget)
        ###***###

        self.comboBox_img = QtWidgets.QComboBox(self.centralwidget)
        self.comboBox_img.setGeometry(QtCore.QRect(690, 590, 401, 31))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        self.comboBox_img.setFont(font)
        self.comboBox_img.setObjectName("comboBox_img")
        self.pushButton_imgR = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_imgR.setGeometry(QtCore.QRect(1110, 590, 51, 31))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(18)
        font.setBold(True)
        font.setWeight(75)
        self.pushButton_imgR.setFont(font)
        self.pushButton_imgR.setObjectName("pushButton_imgR")
        self.pushButton_imgL = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_imgL.setGeometry(QtCore.QRect(620, 590, 51, 31))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(18)
        font.setBold(True)
        font.setWeight(75)
        self.pushButton_imgL.setFont(font)
        self.pushButton_imgL.setObjectName("pushButton_imgL")
        self.label_imagetitle = QtWidgets.QLabel(self.centralwidget)
        self.label_imagetitle.setGeometry(QtCore.QRect(670, 10, 441, 31))
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setUnderline(False)
        self.label_imagetitle.setFont(font)
        self.label_imagetitle.setAlignment(QtCore.Qt.AlignCenter)
        self.label_imagetitle.setObjectName("label_imagetitle")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1269, 21))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        self.menuEdit = QtWidgets.QMenu(self.menubar)
        self.menuEdit.setObjectName("menuEdit")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.actionOpen_CellRanger_Matrix = QtWidgets.QAction(MainWindow)
        self.actionOpen_CellRanger_Matrix.setObjectName("actionOpen_CellRanger_Matrix")
        self.actionOpen_h5ad_session = QtWidgets.QAction(MainWindow)
        self.actionOpen_h5ad_session.setObjectName("actionOpen_h5ad_session")
        self.actiondpi_setting = QtWidgets.QAction(MainWindow)
        self.actiondpi_setting.setObjectName("actiondpi_setting")
        self.actionClear_Console_Output = QtWidgets.QAction(MainWindow)
        self.actionClear_Console_Output.setObjectName("actionClear_Console_Output")
        self.menuFile.addAction(self.actionOpen_CellRanger_Matrix)
        self.menuFile.addAction(self.actionOpen_h5ad_session)
        self.menuEdit.addAction(self.actiondpi_setting)
        self.menuEdit.addAction(self.actionClear_Console_Output)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuEdit.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

        ###***###
        ##finish setting up buttons
        # setting up combobox
        self.comboBox_rankgene.addItem("logreg")
        self.comboBox_rankgene.addItem("wilcoxon")
        # linking buttons to actions
        self.bttn_next.clicked.connect(self.bttn_next_pressed)
        self.bttn_recalc.setText("Calculate")
        self.bttn_recalc.clicked.connect(self.bttn_recalc_pressed)
        self.bttn_excludegene.clicked.connect(self.bttn_excludegene_pressed)
        self.custom_exclusion_text = ''
        self.cycle_path = '-1'
        self.checkBox_cellcyclescore.clicked.connect(self.checkBox_cellcyclescore_clicked)
        self.checkBox_cellcyclescore.setChecked(False)
        self.bttn_cellcyclescore.clicked.connect(self.bttn_cellcyclescore_pressed)
        self.bttn_cellcyclescore.setText("human")
        self.bttn_cellcyclescoretext = "human"
        self.txt_numhivargenes2keep.textChanged.connect(self.txt_numhivargenes2keep_delta)
        self.bttn_NC_listgenes.clicked.connect(self.bttn_NClistgenes_pressed)
        self.NCmarker_genes_txt = ''
        self.cluster_names = ''
        self.bttn_plot.clicked.connect(self.bttn_plot_pressed)
        self.saveplot_on = False
        self.comboBox_img.currentIndexChanged.connect(self.imgChangedComboBox)
        self.pushButton_imgL.pressed.connect(self.imgLaction)
        self.pushButton_imgR.pressed.connect(self.imgRaction)

        # opening file menu to start
        self.actionOpen_CellRanger_Matrix.triggered.connect(self.CellRangerMatrix)
        self.actionOpen_h5ad_session.triggered.connect(self.h5adSession)
        # edit menu
        self.actiondpi_setting.triggered.connect(self.setfigDPI)
        self.dpisetting = 150
        self.whichstage = 0
        self.actionClear_Console_Output.triggered.connect(self.clearConsole)

        # output console (DISABLE THE NEXT 5 LINES IF YOU NEED TO TROUBLESHOOT)
        self.firstcall = True
        XStream.stdout().msgwritten.connect(self.testlogger)
        XStream.stdout().msgwritten.connect(self.plainTextEdit.insertPlainText)
        XStream.stderr().msgwritten.connect(self.testlogger)
        XStream.stderr().msgwritten.connect(self.plainTextEdit.insertPlainText)

        # cache Fig stuff
        self.CacheFigMasterList = ['highest expressed genes (original)',  # 0
                                   'violin: n_gene, n-count, % mito',  # 1
                                   'scatter: n-count vs. % mito',  # 2
                                   'scatter: n-count vs. n-gene',  # 3
                                   'highest expressed genes',  # 4
                                   'pca variance ratio',  # 5
                                   'PCA',  # 6
                                   'filter genes dispersion',  # 7
                                   'umap',  # 8
                                   'rank gene groups - leiden',  # 9
                                   'umap - leiden',  # 10
                                   'umap marker genes',  # 11
                                   'rank gene plot: dot plot',  # 12
                                   'rank gene plot: stacked violin (horiz)',  # 13
                                   'rank gene plot: stacked violin (vert)',  # 14
                                   'rank gene plot: matrix plot',  # 15
                                   'rank gene plot: heatmap',  # 16
                                   'rank gene plot: tracks plot']  # 17
        self.CacheCurIdx = int('0')
        self.CacheFigCurrentList = []


    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "scRNA Basic Tool"))
        self.albl_cell.setText(_translate("MainWindow", "Cells:"))
        self.albl_gene.setText(_translate("MainWindow", "Genes:"))
        self.albl_obsgrp.setText(_translate("MainWindow", "Observe Groups: "))
        self.label.setText(_translate("MainWindow", "Filters:"))
        self.lbl_filtcell.setText(_translate("MainWindow", "Filter cells - expressing min # of genes:"))
        self.lbl_filtgene.setText(_translate("MainWindow", "Filter genes - expression in min # of cells:"))
        self.txt_filtcell.setPlaceholderText(_translate("MainWindow", "200"))
        self.txt_filtgene.setPlaceholderText(_translate("MainWindow", "3"))
        self.bttn_recalc.setText(_translate("MainWindow", "Recalculate"))
        self.bttn_recalc.setShortcut(_translate("MainWindow", "Ctrl+R"))
        self.lbl_genepercell.setText(_translate("MainWindow", "Max genes per cell"))
        self.lbl_totalgenect.setText(_translate("MainWindow", "Max total gene counts"))
        self.lbl_ERCC.setText(_translate("MainWindow", "Max ERCC %"))
        self.lbl_mitoperc.setText(_translate("MainWindow", "Max mito gene %"))
        self.checkBox_stdgenename.setText(_translate("MainWindow", "Standardize Gene Names"))
        self.checkBox_excludegene.setText(_translate("MainWindow", "Exclude mitochondria/ERCC/\n"
" ribosome or other\n"
" genes from analysis"))
        self.checkBox_cellcyclescore.setText(_translate("MainWindow", "Analysis for Cell Cycle Score"))
        self.lbl_numhivargenes2keep.setText(_translate("MainWindow", "# highly variable\n"
" genes to keep"))
        self.txt_numhivargenes2keep.setPlaceholderText(_translate("MainWindow", "0"))
        self.lbl_mmmean_vargenes.setText(_translate("MainWindow", "Min/Max mean for\n"
" highly variable genes"))
        self.txt_minmean_vargenes.setPlaceholderText(_translate("MainWindow", "0.0125"))
        self.txt_maxmean_vargenes.setPlaceholderText(_translate("MainWindow", "5"))
        self.lbl_mmmean_vargenes_slash.setText(_translate("MainWindow", "/"))
        self.lbl_mdispersion_vargenes.setText(_translate("MainWindow", "Min dispersion for\n"
" highly variable genes"))
        self.txt_mdispersion_vargenes.setPlaceholderText(_translate("MainWindow", "0.05"))
        self.checkBox_reg_phase.setText(_translate("MainWindow", "Phase"))
        self.label_2.setText(_translate("MainWindow", "Regress out:"))
        self.checkBox_reg_mito_nct.setText(_translate("MainWindow", "Mito % and n_counts"))
        self.bttn_next.setText(_translate("MainWindow", "Next Param."))
        self.bttn_next.setShortcut(_translate("MainWindow", "Ctrl+N"))
        self.lbl_NC_listgenes.setText(_translate("MainWindow", "List of known genes\n"
" for eval"))
        self.label_3.setText(_translate("MainWindow", "Neighborhood Cluster:"))
        self.bttn_NC_listgenes.setText(_translate("MainWindow", "Type"))
        self.lbl_NC_neighborhoodnum.setText(_translate("MainWindow", "Neighborhood #"))
        self.txt_NC_neighborhoodnum.setPlaceholderText(_translate("MainWindow", "15"))
        self.txt_NC_PCA.setPlaceholderText(_translate("MainWindow", "0"))
        self.lbl_NC_PCA.setText(_translate("MainWindow", "PCA #"))
        self.lbl_NC_resclust.setText(_translate("MainWindow", "Resolution for clustering\n"
" (higher = more clusters)"))
        self.txt_NC_resclust.setPlaceholderText(_translate("MainWindow", "1"))
        self.lbl_NC_dist.setText(_translate("MainWindow", "Min dist."))
        self.txt_NC_dist.setPlaceholderText(_translate("MainWindow", "1"))
        self.txt_NC_dimemb.setPlaceholderText(_translate("MainWindow", "2"))
        self.lbl_NC_dimemb.setText(_translate("MainWindow", "Dim. embedding"))
        self.label_4.setText(_translate("MainWindow", "Rank Gene Groups:"))
        self.bttn_cellcyclescore.setText(_translate("MainWindow", "Type"))
        self.bttn_excludegene.setText(_translate("MainWindow", "Type"))
        self.lbl_rankgene.setText(_translate("MainWindow", "# top ranked genes\n"
" for each group"))
        self.textEdit_genevis.setPlaceholderText(_translate("MainWindow", "default: all (again, use commas)"))
        self.lbl_genevis.setText(_translate("MainWindow", "List genes for visualization:"))
        self.label_5.setText(_translate("MainWindow", "Rank Gene Plotting:"))
        self.lbl_numgeneplot.setText(_translate("MainWindow", "# top ranked genes\n"
" from each group\n"
" to plot:"))
        self.txt_numgeneplot.setPlaceholderText(_translate("MainWindow", "8"))
        self.bttn_plot.setText(_translate("MainWindow", "Generate\n"
"Plot"))
        self.lbl_renameclust.setText(_translate("MainWindow", "Rename clusters:"))
        self.textEdit_renameclust.setPlaceholderText(_translate("MainWindow", "Use commas in between, don\'t hit enter! (otherwise you can leave this blank)"))
        self.checkBox_plot_dot.setText(_translate("MainWindow", "dot plot"))
        self.checkBox_plot_stkviolinH.setText(_translate("MainWindow", "stacked violin (horiz.)"))
        self.checkBox_plot_stkviolinV.setText(_translate("MainWindow", "stacked violin (vert.)"))
        self.checkBox_plot_save.setText(_translate("MainWindow", "Save?"))
        self.checkBox_plot_heat.setText(_translate("MainWindow", "heatmap"))
        self.checkBox_plot_matrix.setText(_translate("MainWindow", "matrix plot"))
        self.checkBox_plot_tracks.setText(_translate("MainWindow", "tracks plot"))
        self.pushButton_imgR.setText(_translate("MainWindow", ">"))
        self.pushButton_imgL.setText(_translate("MainWindow", "<"))
        self.label_imagetitle.setText(_translate("MainWindow", ""))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.menuEdit.setTitle(_translate("MainWindow", "Edit"))
        self.actionOpen_CellRanger_Matrix.setText(_translate("MainWindow", "Open CellRanger Matrix"))
        self.actionOpen_h5ad_session.setText(_translate("MainWindow", "Open h5ad session"))
        self.actiondpi_setting.setText(_translate("MainWindow", "Change Figure DPI"))
        self.actionClear_Console_Output.setText(_translate("MainWindow", "Clear Console Output"))




    ###***###
    #############################################################################################################
    #############################################################################################################
    def CellRangerMatrix(self):
        self.fromCellRangerMatrix=True
        '''
        print(" --- TYPE PROJECT NAME --- ")
        wid1= QInputDialog()
        self.name, okPressed = QInputDialog.getText(wid1, 'Project', 'Enter Project Name:', QLineEdit.Normal, 'Proj1')
        try:
            if okPressed and (self.name != ''):
                print("Project Name: " + self.name)
            else:
                raise ValueError

            wid2= QFileDialog()
            wid2.setFileMode(QFileDialog.Directory)
            curdir=os.getcwd()
            #print(curdir)
            print(" ---  SELECT PROJECT DIRECTORY  --- ")
            self.path = QFileDialog.getExistingDirectory(wid2, 'Select Project Directory', curdir)
            print('proj path: ' + self.path)
            os.chdir(self.path)
            '''
        try:
            self.setProjNamePath()

            curdir = os.getcwd()
            wid3= QFileDialog()
            print(" ---  SELECT DIRECTORY CONTAINING: barcodes.tsv, genes.tsv, matrix.mtx  --- ")
            self.data_path = QFileDialog.getExistingDirectory(wid3, 'Select folder containing barcodes.tsv, genes.tsv, and matrix.mtx', curdir)
            print('data path: ' + self.data_path)
            file_l = [x for x in os.listdir(self.data_path)]
            if  'barcodes.tsv' and 'genes.tsv' and 'matrix.mtx' in file_l:
                print('Correct pathway for matrix')
                # make directory for file saving if not exist
                if not os.path.exists('./write'):
                     os.makedirs('./write')

            else:
                msg = QMessageBox()
                msg.setWindowTitle("Error")
                msg.setText("The barcode, genes, and matrix files are not located in that folder")
                msg.setIcon(QMessageBox.Critical)
                x=msg.exec_()
                print("incorrect directory chosen")
                raise ValueError("incorrect file/directory chosen")

            self.whichstage=0
            self.common_load()
            self.txt_filtcell.setEnabled(True)
            self.txt_filtgene.setEnabled(True)
            self.bttn_next.setEnabled(True)
            self.bttn_recalc.setEnabled(True)
            self.whichstage=1

            self.refreshtext()
            self.refreshComboBoxImg()
            self.label_imagetitle.setText(self.CacheFigCurrentList[self.comboBox_img.currentIndex()])

        except:
            print("something didn't work!\n")

    def h5adSession(self):
        self.fromCellRangerMatrix=False
        try:
            self.setProjNamePath()

            curdir = os.getcwd()
            wid3 = QFileDialog()
            print(" ---  SELECT THE h5ad SESSION  --- ")
            self.data_path, dump = QFileDialog.getOpenFileName(wid3, 'Select the saved h5ad session', curdir, 'h5ad session (*.h5ad)')
            print('data path: ' + self.data_path)
            if '.h5ad' in self.data_path:
                print('>>> Processed from saved session: ', self.data_path.split('/')[-1])
                # make directory for file saving if not exist
                if not os.path.exists('./write'):
                    os.makedirs('./write')
            else:
                msg = QMessageBox()
                msg.setWindowTitle("Error")
                msg.setText("Not a valid h5ad session found")
                msg.setIcon(QMessageBox.Critical)
                x = msg.exec_()
                print("incorrect directory chosen")
                raise ValueError("incorrect file/directory chosen")

            self.whichstage = 0
            self.common_load()
            self.txt_filtcell.setEnabled(True)
            self.txt_filtgene.setEnabled(True)
            self.bttn_next.setEnabled(True)
            self.bttn_recalc.setEnabled(True)
            self.whichstage = 1

            self.refreshtext()
            self.refreshComboBoxImg()
            self.label_imagetitle.setText(self.CacheFigCurrentList[self.comboBox_img.currentIndex()])

        except:
            print("something didn't work!\n")

    def setProjNamePath(self):
        print(" --- TYPE PROJECT NAME --- ")
        wid1 = QInputDialog()
        self.name, okPressed = QInputDialog.getText(wid1, 'Project', 'Enter Project Name:', QLineEdit.Normal, 'Proj1')

        if okPressed and (self.name != ''):
            print("Project Name: " + self.name)
        else:
            raise ValueError

        wid2 = QFileDialog()
        wid2.setFileMode(QFileDialog.Directory)
        curdir = os.getcwd()
        # print(curdir)
        print(" ---  SELECT PROJECT DIRECTORY  --- ")
        self.path = QFileDialog.getExistingDirectory(wid2, 'Select Project Directory', curdir)
        print('proj path: ' + self.path)
        os.chdir(self.path)



    def bttn_recalc_pressed(self):
        self.nxtprs = False
        self.adata = self.adatold.copy()
        try:
            if self.whichstage == 1:
                self.data_filter_p1()
            elif self.whichstage==2:
                self.data_filter_p2()
                self.data_normalize()
            elif self.whichstage==3:
                self.mito_exclusion()
            elif self.whichstage==4:
                self.cell_cycle_score()
                self.PCA()
            elif self.whichstage==5:
                self.highly_variable_genes()
                '''if not(self.cycle_choice.lower() in ['y', 'yes']):
                    self.regress_out()'''
            elif self.whichstage==6:
                self.regress_out()
            elif self.whichstage==7:
                self.neighborhood_cluster()
            elif self.whichstage == 8:
                self.rank_genes_groups()
            elif self.whichstage == 9:
                self.rename_clusters()
                self.adatold = self.adata.copy()
            elif self.whichstage == 10:
                self.rank_genes_plots()

            self.refreshtext()
            self.refreshComboBoxImg()
            self.label_imagetitle.setText(self.CacheFigCurrentList[self.comboBox_img.currentIndex()])
            os.chdir(self.path)
        except:
            os.chdir(self.path)
            print("something went wrong with stage " + str(self.whichstage) + "\n")
            return ValueError("something went wrong with stage" + str(self.whichstage))

    def bttn_next_pressed(self):
        self.nxtprs = True
        self.adata = self.adatold.copy()
        try:
            if self.whichstage==0:
                self.txt_filtcell.setEnabled(True)
                self.txt_filtgene.setEnabled(True)

                self.whichstage = 1
            elif self.whichstage==1:
                self.data_filter_p1()
                self.adatold = self.adata.copy()

                self.txt_filtcell.setEnabled(False)
                self.txt_filtgene.setEnabled(False)
                self.txt_genepercell.setEnabled(True)
                self.txt_totalgenect.setEnabled(True)
                self.txt_mitoperc.setEnabled(True)
                self.txt_ERCC.setEnabled(True)
                self.checkBox_stdgenename.setEnabled(True)

                self.whichstage=2
            elif self.whichstage==2:
                self.data_filter_p2()
                self.data_normalize()
                self.adatold = self.adata.copy()

                self.txt_genepercell.setEnabled(False)
                self.txt_totalgenect.setEnabled(False)
                self.txt_mitoperc.setEnabled(False)
                self.txt_ERCC.setEnabled(False)
                self.checkBox_stdgenename.setEnabled(False)
                self.checkBox_excludegene.setEnabled(True)
                self.bttn_excludegene.setEnabled(True)

                self.whichstage=3
            elif self.whichstage==3:
                self.mito_exclusion()
                self.adatold = self.adata.copy()

                self.checkBox_excludegene.setEnabled(False)
                self.bttn_excludegene.setEnabled(False)
                self.checkBox_cellcyclescore.setEnabled(True)
                self.bttn_cellcyclescore.setEnabled(True)

                self.whichstage=4
            elif self.whichstage==4:
                self.cell_cycle_score()
                self.PCA()
                self.adatold = self.adata.copy()

                self.checkBox_cellcyclescore.setEnabled(False)
                self.bttn_cellcyclescore.setEnabled(False)
                self.txt_numhivargenes2keep.setEnabled(True)
                self.txt_numhivargenes2keep.setPlainText('0')
                self.txt_minmean_vargenes.setEnabled(True)
                self.txt_maxmean_vargenes.setEnabled(True)
                self.txt_mdispersion_vargenes.setEnabled(True)

                self.whichstage=5
            elif self.whichstage==5:
                self.highly_variable_genes()
                self.adatold = self.adata.copy()

                self.txt_numhivargenes2keep.setEnabled(False)
                self.txt_minmean_vargenes.setEnabled(False)
                self.txt_maxmean_vargenes.setEnabled(False)
                self.txt_mdispersion_vargenes.setEnabled(False)
                self.checkBox_reg_mito_nct.setEnabled(True)
                if self.cycle_choice.lower() in ['y', 'yes']:
                    self.checkBox_reg_phase.setEnabled(True)
                else:
                    self.checkBox_reg_phase.setEnabled(False)

                self.whichstage=6
            elif self.whichstage==6:
                self.regress_out()
                self.adatold = self.adata.copy()

                self.checkBox_reg_phase.setEnabled(False)
                self.checkBox_reg_mito_nct.setEnabled(False)
                self.bttn_NC_listgenes.setEnabled(True)
                self.txt_NC_neighborhoodnum.setEnabled(True)
                self.txt_NC_PCA.setEnabled(True)
                self.txt_NC_resclust.setEnabled(True)
                self.txt_NC_dimemb.setEnabled(True)
                self.txt_NC_dist.setEnabled(True)

                self.whichstage=7
            elif self.whichstage==7:
                self.neighborhood_cluster()
                self.adatold = self.adata.copy()

                self.bttn_NC_listgenes.setEnabled(False)
                self.txt_NC_neighborhoodnum.setEnabled(False)
                self.txt_NC_PCA.setEnabled(False)
                self.txt_NC_resclust.setEnabled(False)
                self.txt_NC_dimemb.setEnabled(False)
                self.txt_NC_dist.setEnabled(False)
                self.comboBox_rankgene.setEnabled(True)
                self.txt_rankgenenum.setEnabled(True)

                self.whichstage=8
            elif self.whichstage==8:
                self.rank_genes_groups()
                self.adatold = self.adata.copy()

                self.comboBox_rankgene.setEnabled(False)
                self.txt_rankgenenum.setEnabled(False)
                self.textEdit_genevis.setEnabled(True)
                self.textEdit_renameclust.setEnabled(True)
                self.lbl_renameclust.setText("Rename clusters (" + str(len(self.adata.obs['leiden'].cat.categories)) + "):")

                self.whichstage=9
            elif self.whichstage==9:
                self.rename_clusters()
                self.adatold = self.adata.copy()

                self.textEdit_genevis.setEnabled(False)
                self.textEdit_renameclust.setEnabled(False)
                self.txt_numgeneplot.setEnabled(True)
                self.bttn_plot.setEnabled(True)
                self.checkBox_plot_dot.setEnabled(True)
                self.checkBox_plot_heat.setEnabled(True)
                self.checkBox_plot_stkviolinV.setEnabled(True)
                self.checkBox_plot_stkviolinH.setEnabled(True)
                self.checkBox_plot_matrix.setEnabled(True)
                self.checkBox_plot_tracks.setEnabled(True)
                self.checkBox_plot_save.setEnabled(True)

                self.bttn_next.setEnabled(False)
                self.bttn_recalc.setEnabled(False)

                self.whichstage=10
            elif self.whichstage==10:
                self.rank_genes_plots() #technically should never be able to get here
                #self.adatold = self.adata.copy()

            self.refreshtext()
            self.refreshComboBoxImg()
            self.label_imagetitle.setText(self.CacheFigCurrentList[self.comboBox_img.currentIndex()])
            os.chdir(self.path)
            print("Completed stage " + str(self.whichstage-1) + "\n")
        except:
            os.chdir(self.path)
            print("something went wrong with stage " + str(self.whichstage) + "\n")
            return ValueError



    def bttn_excludegene_pressed(self):
        wid1 = QInputDialog()
        oldtext=self.custom_exclusion_text
        self.custom_exclusion_text, okPressed = QInputDialog.getText(wid1, 'Exclude Genes', 'Enter Genes to Exclude:', QLineEdit.Normal,
                                                                        self.custom_exclusion_text)
        if not okPressed:
            self.custom_exclusion_text=oldtext

    def checkBox_cellcyclescore_clicked(self):
        if self.checkBox_cellcyclescore.isChecked():
            curdir = os.getcwd()
            wid1 = QFileDialog()
            #wid1.setFilter(wid1, "Text files (*.txt)")
            self.cycle_path, test = QFileDialog.getOpenFileName(wid1,'Open file...',curdir,"Text files (*.txt)")
            if self.cycle_path=='':
                self.cycle_path='-1'
                self.checkBox_cellcyclescore.setChecked(False)
            print("Cell Cycle list path: " + self.cycle_path)

    def bttn_cellcyclescore_pressed(self):
        if self.bttn_cellcyclescoretext=="human":
            self.bttn_cellcyclescore.setText("mouse")
            self.bttn_cellcyclescoretext="mouse"
        elif self.bttn_cellcyclescoretext=="mouse":
            self.bttn_cellcyclescore.setText("human")
            self.bttn_cellcyclescoretext="human"

    def txt_numhivargenes2keep_delta(self):
        temp=self.txt_numhivargenes2keep.toPlainText()
        self.varmeanopen=False
        if temp.isdigit():
            temp1=int(temp)
            if temp1<1:
                self.varmeanopen=True
        else:
            self.txt_numhivargenes2keep.setPlainText('0')
            self.varmeanopen=True

        if self.varmeanopen:
            self.txt_minmean_vargenes.setEnabled(True)
            self.txt_maxmean_vargenes.setEnabled(True)
            self.txt_mdispersion_vargenes.setEnabled(True)
        else:
            self.txt_minmean_vargenes.setEnabled(False)
            self.txt_maxmean_vargenes.setEnabled(False)
            self.txt_mdispersion_vargenes.setEnabled(False)

    def bttn_NClistgenes_pressed(self):
        wid1 = QInputDialog()
        oldtext=self.NCmarker_genes_txt
        self.NCmarker_genes_txt, okPressed = QInputDialog.getText(wid1, 'Genes', 'Enter list of known genes to evaluate:', QLineEdit.Normal,
                                                    self.NCmarker_genes_txt)
        if not okPressed:
            self.NCmarker_genes_txt=oldtext

    def bttn_plot_pressed(self):
        self.rank_genes_plots()
        self.refreshComboBoxImg()
        self.label_imagetitle.setText(self.CacheFigCurrentList[self.comboBox_img.currentIndex()])


    ######## Working on Viewer ########
    def refreshViewer(self):
        curboxidx=self.CacheCurIdx
        os.chdir(self.path)
        os.chdir("cache")
        os.chdir("figures")
        temp=os.listdir()
        self.viewer.setImage(QPixmap.fromImage(QImage(temp[curboxidx])))
        #self.label_imagetitle.setText(self.CacheFigCurrentList[self.CacheCurIdx])
        os.chdir(self.path)

    def refreshComboBoxImg(self):
        self.getCacheFigList()
        oldcachefiglist = self.CacheFigCurrentList
        newcachefiglist = []
        newidx=0
        for x in self.CacheFigIdx:
            newcachefiglist.append(self.CacheFigMasterList[x])
        if not (oldcachefiglist == newcachefiglist):
            self.comboBox_img.clear()
            self.CacheFigCurrentList = newcachefiglist
            for y in newcachefiglist:
                self.comboBox_img.addItem(y)

            try:
                self.comboBox_img.setCurrentIndex(len(newcachefiglist)-1)
                newidx = len(newcachefiglist)-1
            except:
                self.comboBox_img.setCurrentIndex(0)
                newidx=0

            self.CacheCurIdx = newidx
            self.refreshViewer()

    def imgChangedComboBox(self):
        self.CacheCurIdx = self.comboBox_img.currentIndex()
        self.refreshViewer()

        templist = self.CacheFigCurrentList
        tempidx = self.CacheCurIdx
        tempname = templist[tempidx]
        self.label_imagetitle.setText(tempname)

    def imgLaction(self):
        curidx=self.CacheCurIdx
        if not(curidx==0) and len(self.CacheFigIdx)>0:
            self.CacheCurIdx=curidx-1
            self.refreshViewer()

            templist = self.CacheFigCurrentList
            tempidx = self.CacheCurIdx
            tempname = templist[tempidx]
            self.label_imagetitle.setText(tempname)
            self.comboBox_img.setCurrentIndex(tempidx)
        else:
            print("you're already on the first figure")

    def imgRaction(self):
        curidx=self.CacheCurIdx
        if not(curidx>=(len(self.CacheFigIdx)-1)) and len(self.CacheFigIdx)>0:
            self.CacheCurIdx=curidx+1
            self.refreshViewer()

            templist = self.CacheFigCurrentList
            tempidx = self.CacheCurIdx
            tempname = templist[tempidx]
            self.label_imagetitle.setText(tempname)
            self.comboBox_img.setCurrentIndex(tempidx)
        else:
            print("you're already on the last figure")

    def getCacheFigList(self):
        os.chdir(self.path)
        os.chdir("cache")
        os.chdir("figures")
        temp = os.listdir()
        tapered = []
        for x in temp:
            tapered.append(int(x[0:2]))
        self.CacheFigIdx=tapered
        os.chdir(self.path)

    ###############################################################################################################
    def refreshtext(self):
        self.albl_cell.setText("Cells: " + str(self.adata.n_obs) + "  (of " + self.cellnum_ori + ")")
        self.albl_gene.setText("Genes: " + str(len(self.adata.var_names)) + "  (of " + self.genenum_ori + ")")
        self.albl_obsgrp.setText("Observe Groups: " + str(len(self.adata.obs.columns.values)))

    def save_data(self):
        print('{} cells x {} genes with {} observe_groups'.format(self.adata.n_obs, len(self.adata.var_names),
                                                                  len(self.adata.obs.columns.values)))
        print(self.adata)
        self.adata.write(self.results_file)
        return 'results_file is saved'

    def parameter_log(self, name, dic):
        parameter_dict = {}
        parameter_dict.update(dic)
        name = str(name)
        with open('./write/parameter_log ' + self.time + '.txt', 'a', encoding='utf-8') as f:
            f.write(name + '\n')
            for k, v in parameter_dict.items():
                f.write(k + ':' + str(v) + '\n')
            f.write('\n')

    def renamereplace(self,strO,strD):
        try:
            os.rename(strO,strD)
        except:
            try:
                os.remove(strD)
                os.rename(strO, strD)
                print("replaced file: "+ strD +" with " +strO)
            except:
                print("why didn't this work?")

    def copyreplace(self,strO,strD):
        try:
            shutil.copyfile(strO,strD)
        except:
            try:
                os.remove(strD)
                shutil.copyfile(strO, strD)
                print("replaced file: "+ strD +" with " +strO)
            except:
                print("why didn't this work?")

    ############################################################################################
    def matrix_load(self):
        self.adata = sc.read_10x_mtx(self.data_path,  # the directory with the `.mtx` file, cell(barcodes.tsv)-row gene(genes.tsv)-column
                                     var_names='gene_symbols',  # use gene symbols for the variable names (variables-axis index)
                                     cache=True)      # write a cache file for faster subsequent reading
        self.adata.var_names_make_unique()  # this is unnecessary if using 'gene_ids'
        #self.adatori = self.adata.copy()
        self.adatold = self.adata.copy()

        self.results_raw_file = './write/{}_raw.h5ad'.format(self.name)  # the file that will store the analysis results
        self.results_file = './write/{}.h5ad'.format(self.name)  # the file that will store the analysis results

        sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
        sc.logging.print_versions()
        sc.settings.set_figure_params(dpi=self.dpisetting, dpi_save=self.dpisetting)
        pd.set_option('display.max_columns', 30)
        self.time=time.strftime("%Y-%m-%d_%H.%M.%S", time.localtime())
        print(self.time)
        print('{} cells x {} genes with {} observe_groups'.format(self.adata.n_obs,len(self.adata.var_names),
                                                                  len(self.adata.obs.columns.values)))
        self.cellnum_ori = str(self.adata.n_obs)
        self.genenum_ori = str(len(self.adata.var_names))
        self.refreshtext()

        # Show those genes with highest (top N) fraction of counts in each single cells, across all cells
        os.chdir("cache")
        try:
            shutil.rmtree("figures")
        except:
            dump = 1
        try:
            os.mkdir("figures")
        except:
            print("make sure you're not in the ./cache/figures folder please!")
            raise ValueError
        sc.pl.highest_expr_genes(self.adata, n_top=20, save='_{}_{}.png'.format(self.time, self.name))
        self.renamereplace('./figures/highest_expr_genes_' + self.time + '_' + self.name + '.png',
                  './figures/00_' + self.name + '_' + self.time + '_highest_expr_genes_original.png')
        # os.rename('highest_expr_genes_.png','001_highest_expr_genes.png')
        os.chdir("..")

    def session_load(self):
        self.adata = sc.read(self.data_path)
        self.adata.var_names_make_unique()  # this is unnecessary if using 'gene_ids'
        #self.adatori = self.adata.copy()
        self.adatold = self.adata.copy()

        self.results_raw_file = './write/{}_raw.h5ad'.format(self.name)  # the file that will store the analysis results
        self.results_file = './write/{}.h5ad'.format(self.name)  # the file that will store the analysis results


    def common_load(self):
        if self.fromCellRangerMatrix:
            self.adata = sc.read_10x_mtx(self.data_path, # the directory with the `.mtx` file, cell(barcodes.tsv)-row gene(genes.tsv)-column
                                         var_names='gene_symbols', # use gene symbols for the variable names (variables-axis index)
                                         cache=True)  # write a cache file for faster subsequent reading
        else:
            self.adata = sc.read(self.data_path)
            os.chdir(self.path)
            try:
                os.chdir("cache")
                os.chdir("..")
            except:
                os.mkdir("cache")
        self.adata.var_names_make_unique()  # this is unnecessary if using 'gene_ids'

        self.results_raw_file = './write/{}_raw.h5ad'.format(self.name)  # the file that will store the analysis results
        self.results_file = './write/{}.h5ad'.format(self.name)  # the file that will store the analysis results
        sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
        sc.logging.print_versions()
        sc.settings.set_figure_params(dpi=self.dpisetting)
        pd.set_option('display.max_columns', 30)
        self.time = time.strftime("%Y-%m-%d_%H.%M.%S", time.localtime())
        print(self.time)
        print('{} cells x {} genes with {} observe_groups'.format(self.adata.n_obs, len(self.adata.var_names),
                                                                  len(self.adata.obs.columns.values)))
        self.cellnum_ori = str(self.adata.n_obs)
        self.genenum_ori = str(len(self.adata.var_names))
        self.refreshtext()

        # Show those genes with highest (top N) fraction of counts in each single cells, across all cells
        os.chdir("cache")
        try:
            shutil.rmtree("figures")
        except:
            dump = 1
        try:
            os.mkdir("figures")
        except:
            print("make sure you're not in the ./cache/figures folder please!")
            raise ValueError
        sc.pl.highest_expr_genes(self.adata, n_top=20, save='_{}_{}.png'.format(self.time, self.name))
        self.renamereplace('./figures/highest_expr_genes_' + self.time + '_' + self.name + '.png',
                           './figures/00_' + self.name + '_' + self.time + '_highest_expr_genes_original.png')
        # os.rename('highest_expr_genes_.png','001_highest_expr_genes.png')
        os.chdir("..")

        #self.adatori = self.adata.copy()
        self.adatold = self.adata.copy()



    ######################
    def data_filter_p1(self):

        print('Standard filter')
        try:
            if not(self.txt_filtcell.toPlainText().isdigit() and self.txt_filtgene.toPlainText().isdigit()):
                self.txt_filtcell.setPlainText("200")
                self.txt_filtgene.setPlainText("3")

            self.min_genes = int(self.txt_filtcell.toPlainText())
            self.min_cells = int(self.txt_filtgene.toPlainText())

            sc.pp.filter_cells(self.adata, min_genes=self.min_genes)
            sc.pp.filter_genes(self.adata, min_cells=self.min_cells)
            #self.adata=sc.pp.filter_cells(self.adata, min_genes=self.min_genes, copy=True) #copy=True allows to recalculate? but may be depreciated in future version
            #self.adata=sc.pp.filter_genes(self.adata, min_cells=self.min_cells, copy=True) #copy=True allows to recalculate? but may be depreciated in future version
            print('{} cells x {} genes with {} observe_groups'.format(self.adata.n_obs, len(self.adata.var_names),
                                                                      len(self.adata.obs.columns.values)))
            self.refreshtext()

            sc.pp.calculate_qc_metrics
            # adata from csv. is different from mtx. (use dataframe analysis)
            self.mito_genes = [x for x in self.adata.var_names if (x[0:2] == 'MT' or 'chrM' in x)]
            self.ercc = [x for x in self.adata.var_names if 'ERCC' in x]
            self.ribo_genes = [x for x in self.adata.var_names if ((x[0:3].upper()=='RPL') or (x[0:3].upper()=='RPS'))]
            # for each cell (.X) compute fraction of counts in mito genes vs. all genes and save in the array
            # the `.A1` or np.ravel() is only necessary as X is sparse (to transform to a dense array after summing)
            # add the obs (observation) data as row into adata
            self.adata.obs['percent_mito'] = np.ravel(np.sum(self.adata[:, self.mito_genes].X, axis=1)) / np.ravel(np.sum(self.adata.X, axis=1))
            self.adata.obs['percent_ercc'] = np.ravel(np.sum(self.adata[:, self.ercc].X, axis=1)) / np.ravel(np.sum(self.adata.X, axis=1))
            # add the total counts per cell as observations-annotation to adata
            self.adata.obs['n_counts'] = np.ravel(self.adata.X.sum(axis=1))

            #plotting
            os.chdir("cache")
            sc.pl.violin(self.adata, ['n_genes', 'n_counts', 'percent_mito'], jitter=0.4, multi_panel=True, save='_{}_{}.png'.format(self.time,self.name))
            self.renamereplace('./figures/violin_'+self.time+'_'+self.name+'.png','./figures/01_'+self.name+'_'+self.time+'_violin1.png')
            # plot the data and filter unwanted cells based on plot
            sc.pl.scatter(self.adata, x='n_counts', y='percent_mito', save='_{}_{}.png'.format(self.time,self.name))
            self.renamereplace('./figures/scatter_'+self.time+'_'+self.name+'.png','./figures/02_'+self.name+'_'+self.time+'_scatter1.png')
            # sc.pl.scatter(self.adata, x='n_counts', y='percent_ercc')
            sc.pl.scatter(self.adata, x='n_counts', y='n_genes', save='_{}_{}.png'.format(self.time,self.name))
            self.renamereplace('./figures/scatter_'+self.time+'_'+self.name+'.png','./figures/03'+self.name+'_'+self.time+'_scatter2.png')
            os.chdir("..")

        except:
            print("something didn't work!\n")
            raise ValueError


    ######################
    def data_filter_p2(self):
        try:
            if self.txt_genepercell.toPlainText().isdigit and self.txt_totalgenect.toPlainText().isdigit:
                self.max_n_genes = int(self.txt_genepercell.toPlainText())
                self.max_n_counts = int(self.txt_totalgenect.toPlainText())
            else:
                raise ValueError("Max genes per cell and/or Max total gene counts - not positive integers!")

            tempmito = float(self.txt_mitoperc.toPlainText())
            tempercc = float(self.txt_ERCC.toPlainText())
            if tempmito < 0 or tempercc < 0:
                raise ValueError("Mito gene or ERCC percentages must be 0 or above")
            else:
                self.max_percent_mito = tempmito
                self.max_percent_ercc = tempercc

            self.adata = self.adata[self.adata.obs['n_genes'] < self.max_n_genes, :]
            self.adata = self.adata[self.adata.obs['percent_mito'] < self.max_percent_mito, :]
            self.adata = self.adata[self.adata.obs['percent_ercc'] < self.max_percent_ercc, :]
            self.adata = self.adata[self.adata.obs['n_counts'] < self.max_n_counts, :]

            print('{} cells x {} genes with {} observe_groups'.format(self.adata.n_obs, len(self.adata.var_names),
                                                                  len(self.adata.obs.columns.values)))
            self.refreshtext()

            if self.checkBox_stdgenename.isChecked():
                # Get gene id without Chr
                gene_id = [y.split('_')[0] for y in self.adata.var_names]
                self.mito_genes = [y.split('_')[0] for y in self.mito_genes]
                self.ercc = [y.split('_')[0] for y in self.ercc]
                self.ribo_genes = [y.split('_')[0] for y in self.ribo_genes]
                self.adata.var_names = gene_id
                os.chdir("cache")
                sc.pl.highest_expr_genes(self.adata, n_top=20, save='_{}_{}.png'.format(self.time,self.name))
                self.renamereplace('./figures/highest_expr_genes_'+self.time+'_'+self.name+'.png','./figures/04_'+self.name+'_'+self.time+'_highest_expr_genes2.png')
                os.chdir("..")

            if self.nxtprs:
                self.data_filter_parameter = {'min_genes': self.min_genes, 'min_cells': self.min_cells,
                                              'max_n_genes': self.max_n_genes,'max_n_counts': self.max_n_counts,
                                              'max_percent_mito': self.max_percent_mito}
                self.parameter_log('data_filter_parameter', self.data_filter_parameter)
                self.save_data()

        except:
            print("something didn't work!\n")
            raise ValueError

    def data_normalize(self):
        # default normalization by 10,000 and log-transforms the result
        sc.pp.normalize_per_cell(self.adata, counts_per_cell_after=1e4) #copy=True allows to recalculate? but may be depreciated in future version)
        sc.pp.log1p(self.adata)
        #self.adata=sc.pp.log1p(self.adata) #copy=True allows to recalculate? but may be depreciated in future version
        # Set the .raw attribute to the logarithmized raw gene expression for later use in differential testing and visualizations
        if self.nxtprs:
            self.adata.raw = self.adata  # freezes state of AnnData object - normalized, logarithmized, not corrected
            self.adata.write(self.results_raw_file)
            self.save_data()

    ######################
    def mito_exclusion(self):
        try:
            os.chdir("cache")
            sc.pl.highest_expr_genes(self.adata, n_top=20, save='_{}_{}.png'.format(self.time,self.name))
            self.renamereplace('./figures/highest_expr_genes_'+self.time+'_'+self.name+'.png','./figures/04_'+self.name+'_'+self.time+'_highest_expr_genes2.png')
            os.chdir("..")

            if self.checkBox_excludegene.isChecked():
                self.custom_exclusion = [x for x in self.adata.var_names if x in re.split(',|_|-|!| ', self.custom_exclusion_text)]
                gene_inclusion = [x for x in self.adata.var_names if x not in (self.mito_genes + self.ercc + self.ribo_genes + self.custom_exclusion)]
                self.adata = self.adata[:, gene_inclusion]
                exclusion = [x for x in self.adata.var_names if x in (self.mito_genes + self.ercc + self.ribo_genes + self.custom_exclusion)]
                print('Mitochondrial/ERCC genes are excluded'+'\n'+'{} from input are excluded'.format(exclusion))
                os.chdir("cache")
                sc.pl.highest_expr_genes(self.adata, n_top=20, save='_{}_{}.png'.format(self.time,self.name))
                self.renamereplace('./figures/highest_expr_genes_'+self.time+'_'+self.name+'.png','./figures/04_'+self.name+'_'+self.time+'_highest_expr_genes2.png')
                os.chdir("..")
            else:
                exclusion = ''
                print('Mitochondrial/ERCC genes are included')
            self.refreshtext()
            if self.nxtprs:
                gene_exclusion = {'exclusion': exclusion}
                self.parameter_log('gene_exclusion', gene_exclusion)
                self.save_data()
        except:
            print("something didn't work!\n")
            raise ValueError

    ######################
    def cell_cycle_score(self):
        #self.cycle_choice = input('Analysis for cell cycle score? You need cycle gene list in text file. y for yes, any other key for no\n')
        if self.cycle_path=='-1':
            self.cycle_choice='n'
            cell_cycle_dict = 'no cell cycle score'
            return 'No cell cycle computed'
        else:
            self.cycle_choice='y'
            cycle_path = self.cycle_path #input('Provide file for cell cycle genes:\n').strip("'")
            source = self.bttn_cellcyclescoretext
            if source in ['m', 'mouse']:
                source = 'mouse'
                cell_cycle_genel = [(x.strip()[0] + x.strip()[1:].lower()) for x in open(cycle_path)]
                cell_cycle_genes = [x for x in cell_cycle_genel if x in self.adata.var_names]
            else:
                source = 'human'
                cell_cycle_genel = [x.strip()[0] + x.strip()[1:] for x in open(cycle_path)]
                cell_cycle_genes = [x for x in cell_cycle_genel if x in self.adata.var_names]
            s_genes = [x[1] for x in enumerate(cell_cycle_genes) if cell_cycle_genel.index(x[1]) <= 44]
            g2m_genes = [x[1] for x in enumerate(cell_cycle_genes) if cell_cycle_genel.index(x[1]) > 44]

            cell_cycle_dict = {'s_genes': s_genes, 'g2m_genes': g2m_genes}
            sc.tl.score_genes_cell_cycle(self.adata, s_genes=s_genes, g2m_genes=g2m_genes)


        self.refreshtext()

        if self.nxtprs:
            self.parameter_log('cell_cycle_score', cell_cycle_dict)
            self.save_data()

    ######################
    def PCA(self):
        # PCA for cell cycle before high variable selection
        sc.tl.pca(self.adata, svd_solver='arpack')
        os.chdir("cache")
        sc.pl.pca_variance_ratio(self.adata, log=True, save='_{}_{}.png'.format(self.time,self.name))
        self.renamereplace('./figures/pca_variance_ratio_'+self.time+'_'+self.name+'.png','./figures/05_'+self.name+'_'+self.time+'_pca_variance_ratio.png')
        os.chdir("..")
        if self.nxtprs:
            self.save_data()

    ######################
    def highly_variable_genes(self):
        # Identify highly-variable genes
        #print('To enrich highly variable genes for analysis, provide following parameters \nYou can adjust parameters based on returned graph')
        try:
            #n_top_genes = int(input('Number of highly-variable genes to keep ("0" for setting other parameters): '))
            n_top_genes = int(self.txt_numhivargenes2keep.toPlainText())

            if n_top_genes != 0:
                sc.pp.highly_variable_genes(self.adata, n_top_genes=n_top_genes)
                min_mean, max_mean, min_disp = None, None, None
            else:
                #min_mean = float(input('Min mean for highly variable genes(try 0.0125): '))
                #max_mean = float(input('Max mean for highly variable genes(try 5): '))
                #min_disp = float(input('Min disperation for highly variable genes(try 0.05): '))

                if self.txt_minmean_vargenes.toPlainText()=='':
                    self.txt_minmean_vargenes.setPlainText("0.0125")
                if self.txt_maxmean_vargenes.toPlainText()=='':
                    self.txt_maxmean_vargenes.setPlainText("5")
                if self.txt_mdispersion_vargenes.toPlainText()=='':
                    self.txt_mdispersion_vargenes.setPlainText("0.05")

                try:
                    min_mean=float(self.txt_minmean_vargenes.toPlainText())
                except:
                    self.txt_minmean_vargenes.setPlainText("0.0125")
                    min_mean = float(self.txt_minmean_vargenes.toPlainText())

                try:
                    max_mean=float(self.txt_maxmean_vargenes.toPlainText())
                except:
                    self.txt_maxmean_vargenes.setPlainText("5")
                    max_mean = float(self.txt_maxmean_vargenes.toPlainText())

                try:
                    min_disp=float(self.txt_mdispersion_vargenes.toPlainText())
                except:
                    self.txt_mdispersion_vargenes.setPlainText("0.05")
                    min_disp = float(self.txt_mdispersion_vargenes.toPlainText())
                    #print("your variable parameters have to be floats!")

                if not(min_mean>=0 and max_mean>=0 and min_disp>=0):
                    ValueError('your variable parameters have to be positive!')

                sc.pp.highly_variable_genes(self.adata, min_mean=min_mean, max_mean=max_mean, min_disp=min_disp)
                n_top_genes = None

            self.var_min_mean = min_mean
            self.var_max_mean = max_mean
            self.var_min_disp = min_disp
            #if not(self.cycle_choice.lower() in ['y', 'yes']):
            os.chdir("cache")
            sc.pl.highly_variable_genes(self.adata, save='_{}_{}.png'.format(self.time,self.name))
            self.renamereplace('./figures/filter_genes_dispersion_'+self.time+'_'+self.name+'.png','./figures/07_'+self.name+'_'+self.time+'_filter_genes_dispersion.png')
            os.chdir("..")

            adata_variable = self.adata[:, self.adata.var['highly_variable']]
            print('Parameters for this plot: n_top_genes={},min_mean={},max_mean={},min_disp={}'.format(n_top_genes, min_mean, max_mean, min_disp))
            print('{} cells x {} genes with {} observe_groups'.format(adata_variable.n_obs, len(adata_variable.var_names), len(adata_variable.obs.columns.values)))
            self.refreshtext()

            #if self.cycle_choice.lower() in ['y', 'yes']:
            #copy over from regress_out to get PCA plot
            os.chdir("cache")
            sc.pl.pca(self.adata, color='phase', save='_{}_{}.png'.format(self.time, self.name))
            self.renamereplace('./figures/pca_' + self.time + '_' + self.name + '.png', './figures/06_' + self.name + '_' + self.time + '_pca.png')
            os.chdir("..")
        except:
            return ValueError('something went wrong!\n')

        if self.nxtprs:
            self.adata=self.adata[:, self.adata.var['highly_variable']].copy()
            self.save_data()
            self.refreshtext()
    ######################
    def regress_out(self):
        #choice = input('Based on phase embedded in PCA plot, do you want to regress out phase? y for yes, any other key for no\n')
        #choice_2 = input('Do you want to regree out mito_percentage and n_counts? y for yes, any other key for no\n')
        if self.checkBox_reg_phase.isChecked() and self.checkBox_reg_mito_nct.isChecked():
            regress_out_comp = ['n_counts', 'percent_mito', 'S_score', 'G2M_score']
            sc.pp.regress_out(self.adata, regress_out_comp)
        elif not(self.checkBox_reg_phase.isChecked()) and self.checkBox_reg_mito_nct.isChecked():
            regress_out_comp = ['n_counts', 'percent_mito']
            sc.pp.regress_out(self.adata, regress_out_comp)
        elif self.checkBox_reg_phase.isChecked() and not(self.checkBox_reg_mito_nct.isChecked()):
            regress_out_comp = ['S_score', 'G2M_score']
            sc.pp.regress_out(self.adata, regress_out_comp)
        else:
            regress_out_comp=''
            print('No regression out has been conducted!')

        self.refreshtext()
        self.refreshComboBoxImg()
        self.label_imagetitle.setText(self.CacheFigCurrentList[self.comboBox_img.currentIndex()])

        if self.nxtprs:
            while True:
                wid1= QInputDialog()
                tmaxvalue, okPressed = QInputDialog.getText(wid1, 'Set Scale', 'Based on the graph, set the max value for scale:', QLineEdit.Normal, '10')
                try:
                    if okPressed and (tmaxvalue != ''):
                        max_value = float(tmaxvalue) #float(input('Based on the graph, set the max value for scale: '))
                        break
                    else:
                        raise ValueError("set the value of the scale! Come on!")
                except:
                    print('try again!')

            sc.pp.scale(self.adata, max_value=max_value) #copy=True allows to recalculate? but may be depreciated in future version)
            self.highly_variable_parameter = {'min_mean': self.var_min_mean, 'max_mean': self.var_max_mean, 'min_disp': self.var_min_disp, 'max_value': max_value, 'regress_out': regress_out_comp}

            self.parameter_log('highly_variable/regression_parameter', self.highly_variable_parameter)
            self.save_data()
            self.refreshtext()

    ######################
    def neighborhood_cluster(self):
        print('\n\nCompute neighborhood graph of cells using the PCA representation of the data matrix')
        print('Embed neighborhood graph of cells using UMAP')
        print('Clustering of cells using leiden')
        print('To evaluate embedding and clustering, prepare known marker genes for expected clusters')


        #marker_gene = input('List of known genes for evaluation: \n')
        marker_gene = [x for x in self.adata.raw.var_names if x in re.split(',|_|-|!| ', self.NCmarker_genes_txt)]
        print("\nComputing neighborhood graph")
        print(marker_gene)
        if not(self.txt_NC_neighborhoodnum.toPlainText().isdigit()):
            self.txt_NC_neighborhoodnum.setPlainText("15")
        if not(self.txt_NC_PCA.toPlainText().isdigit()):
            self.txt_NC_PCA.setPlainText("0")
        if not(self.txt_NC_dimemb.toPlainText().isdigit()):
            self.txt_NC_dimemb.setPlainText("2")
        if not(self.txt_NC_dist.toPlainText().isdigit()):
            self.txt_NC_dist.setPlainText("1")
        n_neighbors = int(self.txt_NC_neighborhoodnum.toPlainText()) #int(input('Neighbor numbers for computation (default 15): '))
        n_pcs = int(self.txt_NC_PCA.toPlainText()) #int(input('PCA numbers for computation (default None, input 0): '))
        try:
            leiden_resolution = float(self.txt_NC_resclust.toPlainText())
        except:
            self.txt_NC_resclust.setPlainText("1")
            leiden_resolution = float(self.txt_NC_resclust.toPlainText())
        #leiden_resolution = float(input('Resolution for clustering, higher with more and smaller clusters (default 1): '))
        umap_n_components = int(self.txt_NC_dimemb.toPlainText())#int(input('Dimensions of the embedding (default 2): '))
        try:
            umap_min_dist = float(self.txt_NC_dist.toPlainText())
        except:
            self.txt_NC_dist.setPlainText("1")
            umap_min_dist = float(self.txt_NC_resclust.toPlainText())
        #umap_min_dist = float(input('Effective minimum distance between embedded points, larger values will result on a more even dispersal of points (default 1): '))
        if n_pcs == 0:
            n_pcs = None

        sc.pp.neighbors(self.adata, n_neighbors=n_neighbors, n_pcs=n_pcs)  # neighborhood parameter
        sc.tl.leiden(self.adata)#, resolution=leiden_resolution)
        sc.tl.umap(self.adata, n_components=umap_n_components, min_dist=umap_min_dist)
        os.chdir("cache")
        sc.pl.umap(self.adata, color=['leiden'] + marker_gene, projection='2d', ncols=5, legend_loc='on data', palette='Set1', legend_fontsize=15, size=20, save='_{}_{}.png'.format(self.time,self.name))
        self.renamereplace('./figures/umap_'+self.time+'_'+self.name+'.png','./figures/08_'+self.name+'_'+self.time+'_umap.png')
        os.chdir("..")

        print('Parameters for this plot:\n')
        print('n_neighbors={},n_pcs={},leiden_resolution={},umap_n_components={},umap_min_dist={},'.format(n_neighbors, n_pcs, leiden_resolution, umap_n_components, umap_min_dist))
        self.refreshtext()
        if self.nxtprs:
            self.cluster_parameter = {'n_neighbors': n_neighbors, 'n_pcs': n_pcs, 'leiden_resolution': leiden_resolution,
                                      'umap_n_components': umap_n_components, 'umap_min_dist': umap_min_dist}
            self.parameter_log('cluster_parameter', self.cluster_parameter)
            self.save_data()

    ######################
    def rank_genes_groups(self):
        # find rank_genes_groups and marker genes for clustered group (logreg)
        # logreg give best cluster identification *only scores but no p-values or logfoldchange, consider Diffxpy
        method1 = self.comboBox_rankgene.currentText()

        sc.tl.rank_genes_groups(self.adata, 'leiden', method=method1)  # method:t-test,wilcoxon,logreg(no pvals or fc)
        os.chdir("cache")
        sc.pl.rank_genes_groups(self.adata, n_genes=25, sharey=False, save='_{}_{}.png'.format(self.time,self.name))
        self.renamereplace('./figures/rank_genes_groups_leiden_'+self.time+'_'+self.name+'.png','./figures/09_'+self.name+'_'+self.time+'_rank_genes_groups_leiden.png')
        os.chdir("..")

        self.adata.write(self.results_file)  # write into result
        m = pd.DataFrame(self.adata.uns['rank_genes_groups']['names']).head(50)
        print(m)

        self.marker_genes = list()
        self.marker_gene_dict = {}

        if not(self.txt_rankgenenum.toPlainText().isdigit()):
            self.txt_rankgenenum.setPlainText("25")
        n_genes = int(self.txt_rankgenenum.toPlainText()) #int(input('Number of top ranked genes for each group: '))


        for i in range(m.shape[1]):
            marker_gene_l = []
            for j in range(n_genes):
                marker_gene_l.append(m.iloc[j, i])
                if m.iloc[j, i] not in self.marker_genes:
                    # Add column i, row j element into list if it's unique, .iloc[,] returns series while iloc[,][] returns value
                    self.marker_genes.append(m.iloc[j, i])
            self.marker_gene_dict[m.columns.tolist()[i]] = marker_gene_l
        print(self.marker_gene_dict)
        self.refreshtext()

        if self.nxtprs:
            result = self.adata.uns['rank_genes_groups']
            groups = result['names'].dtype.names
            pd.set_option('display.max_columns', 30)
            pd.set_option('precision', 3)
            df = pd.DataFrame(
                {group + '_' + key[:1]: result[key][group]
                 for group in groups for key in ['names', 'scores']}).head(1000)
            with pd.ExcelWriter('{}_rank_gene_groups.xlsx'.format(self.name)) as writer:
                df.to_excel(writer, index=None, sheet_name='method')
            writer.save()
            writer.close()

            self.rank_genes_parameter = {'Comparison method': method1,
                                         'Number of top ranked genes for each group': n_genes}
            self.parameter_log('rank_genes_parameter', self.rank_genes_parameter)
            self.save_data()

    ######################
    def rename_clusters(self):
        print('\nTop ranked genes for each cluster:\n', self.marker_gene_dict)
        marker_gene = self.textEdit_genevis.toPlainText() #input('List of genes for visualization: \n')
        #marker_gene = [x for x in self.adata.var_names if x in re.split(',|_|-|!| ', marker_gene)]
        if not(marker_gene.isalnum()) or marker_gene=="all":
            self.textEdit_genevis.setPlainText("all")
            marker_gene = self.textEdit_genevis.toPlainText()
            sc.pl.umap(self.adata, color=['leiden'] + self.marker_genes, projection='2d',
                       ncols=5, legend_loc='on data', palette='Set1', legend_fontsize=15, size=20, save='_leiden.png')
        else:
            marker_gene = [x for x in self.adata.var_names if x in re.split(',|_|-|!| ', marker_gene)]
            sc.pl.umap(self.adata, color=['leiden'] + marker_gene, projection='2d',
                       ncols=5, legend_loc='on data', palette='Set1', legend_fontsize=15, size=20, save='_leiden.png')
        self.copyreplace('./figures/umap_leiden.png','./cache/figures/10_'+self.name+'_'+self.time+'_umap_leiden.png')
        try:
            cluster_name = self.textEdit_renameclust.toPlainText() #input('Rename the clusters based on top ranked genes: ')
            cluster_name = re.split(',|_|-|!| ', cluster_name)
            self.adata.rename_categories('leiden', cluster_name)
        except:
            print(('Cluster numbers should be ' + str(len(self.adata.obs['leiden'].cat.categories))))
            return ValueError
        #marker_gene = input('List of genes for visualization, "all" for all marker genes: \n')
        if marker_gene.lower() != 'all':
            marker_gene = [x for x in self.adata.raw.var_names if x in re.split(',|_|-|!| ', marker_gene)]
            sc.pl.umap(self.adata, color=['leiden'] + marker_gene, projection='2d', palette='Set1',
                       legend_fontsize=15, size=20, save='_selected_genes.png')
        else:
            sc.pl.umap(self.adata, color=['leiden'] + self.marker_genes, projection='2d', palette='Set1',
                       legend_fontsize=15, size=20, save='_marker_genes.png')
        self.copyreplace('./figures/umap_marker_genes.png','./cache/figures/11_' + self.name + '_' + self.time + '_umap_marker_genes.png')
        if self.nxtprs:
            self.save_data()

    ######################
    def rank_genes_plots(self):

        if not(self.txt_numgeneplot.toPlainText().isdigit()):
            self.txt_numgeneplot.setPlainText("8")

        n_genes = int(self.txt_numgeneplot.toPlainText()) #int(input('Number of top ranked genes from each group for plot: '))
        self.adata = sc.read(self.results_file)
        if self.checkBox_plot_save.isChecked():
            print("Saving rank_gene_plot figures...")
            # Visualize marker genes in clusters
            if self.checkBox_plot_dot.isChecked():
                sc.pl.rank_genes_groups_dotplot(self.adata, n_genes=n_genes, save='_{}_rank_genes_group.png'.format(self.name))
                self.copyreplace('./figures/dotplot_{}_rank_genes_group.png'.format(self.name),
                                './cache/figures/12_' + self.name + '_' + self.time + '_rank_genes_group_dotplot.png')
            if self.checkBox_plot_stkviolinH.isChecked():
                sc.pl.rank_genes_groups_stacked_violin(self.adata, n_genes=n_genes, save='_{}_rank_genes_group.png'.format(self.name))
                self.copyreplace('./figures/stacked_violin_{}_rank_genes_group.png'.format(self.name),
                                './cache/figures/13_' + self.name + '_' + self.time + '_rank_genes_group_stacked_violinH.png')
            if self.checkBox_plot_stkviolinV.isChecked():
                sc.pl.rank_genes_groups_stacked_violin(self.adata, n_genes=n_genes, swap_axes=True, save='_{}_rank_genes_group-swap.png'.format(self.name))
                self.copyreplace('./figures/stacked_violin_{}_rank_genes_group-swap.png'.format(self.name),
                                './cache/figures/14_' + self.name + '_' + self.time + '_rank_genes_group_stacked_violinV.png')
            if self.checkBox_plot_matrix.isChecked():
                sc.pl.rank_genes_groups_matrixplot(self.adata, n_genes=n_genes, cmap='bwr', save='_{}_rank_genes_group.png'.format(self.name))
                self.copyreplace('./figures/matrixplot_{}_rank_genes_group.png'.format(self.name),
                                './cache/figures/15_' + self.name + '_' + self.time + '_rank_genes_group_matrixplot.png')
            if self.checkBox_plot_heat.isChecked():
                sc.pl.rank_genes_groups_heatmap(self.adata, n_genes=n_genes, swap_axes=True, save='_{}_rank_genes_group.png'.format(self.name))
                self.copyreplace('./figures/heatmap_{}_rank_genes_group.png'.format(self.name),
                                './cache/figures/16_' + self.name + '_' + self.time + '_rank_genes_group_heatmap.png')

            # visualize the tracksplot
            ad = self.adata.copy()  # track plot data is better visualized using the non-log counts
            ad.raw.X.data = np.expm1(ad.raw.X.data)
            if self.checkBox_plot_tracks.isChecked():
                sc.pl.rank_genes_groups_tracksplot(ad, n_genes=n_genes, save='_{}_rank_genes_group.png'.format(self.name))
                self.copyreplace('./figures/tracksplot_{}_rank_genes_group.png'.format(self.name),
                                 './cache/figures/17_' + self.name + '_' + self.time + '_rank_genes_group_tracksplot.png')
        else:
            # Visualize marker genes in clusters
            if self.checkBox_plot_dot.isChecked():
                os.chdir("cache")
                sc.pl.rank_genes_groups_dotplot(self.adata, n_genes=n_genes, save='_{}_{}.png'.format(self.time,self.name))
                self.renamereplace('./figures/dotplot_'+self.time+'_'+self.name+'.png','./figures/12_' + self.name + '_' + self.time + '_rank_genes_group_dotplot.png')#, save='_{}_rank_genes_group.png'.format(self.name))
                os.chdir("..")
            if self.checkBox_plot_stkviolinH.isChecked():
                os.chdir("cache")
                sc.pl.rank_genes_groups_stacked_violin(self.adata, n_genes=n_genes, save='_{}_{}.png'.format(self.time,self.name))
                self.renamereplace('./figures/stacked_violin_'+self.time+'_'+self.name+'.png','./figures/13_' + self.name + '_' + self.time + '_rank_genes_group_stacked_violinH.png')#, save='_{}_rank_genes_group.png'.format(self.name))
                os.chdir("..")
            if self.checkBox_plot_stkviolinV.isChecked():
                os.chdir("cache")
                sc.pl.rank_genes_groups_stacked_violin(self.adata, n_genes=n_genes, swap_axes=True, save='_{}_{}.png'.format(self.time,self.name))
                self.renamereplace('./figures/stacked_violin_'+self.time+'_'+self.name+'.png','./figures/14_' + self.name + '_' + self.time + '_rank_genes_group_stacked_violinV.png')#, save='_{}_rank_genes_group-swap.png'.format(self.name))
                os.chdir("..")
            if self.checkBox_plot_matrix.isChecked():
                os.chdir("cache")
                sc.pl.rank_genes_groups_matrixplot(self.adata, n_genes=n_genes, cmap='bwr', save='_{}_{}.png'.format(self.time,self.name))
                self.renamereplace('./figures/matrixplot_'+self.time+'_'+self.name+'.png','./figures/15_' + self.name + '_' + self.time + '_rank_genes_group_matrixplot.png')#, save='_{}_rank_genes_group.png'.format(self.name))
                os.chdir("..")
            if self.checkBox_plot_heat.isChecked():
                os.chdir("cache")
                sc.pl.rank_genes_groups_heatmap(self.adata, n_genes=n_genes, swap_axes=True, save='_{}_{}.png'.format(self.time,self.name))
                self.renamereplace('./figures/heatmap_'+self.time+'_'+self.name+'.png','./figures/16_' + self.name + '_' + self.time + '_rank_genes_group_heatmap.png')#, save='_{}_rank_genes_group.png'.format(self.name))
                os.chdir("..")

            # visualize the tracksplot
            ad = self.adata.copy()  # track plot data is better visualized using the non-log counts
            ad.raw.X.data = np.expm1(ad.raw.X.data)
            if self.checkBox_plot_tracks.isChecked():
                os.chdir("cache")
                sc.pl.rank_genes_groups_tracksplot(ad, n_genes=n_genes, save='_{}_{}.png'.format(self.time,self.name))
                self.renamereplace('./figures/tracksplot_'+self.time+'_'+self.name+'.png','./figures/17_' + self.name + '_' + self.time + '_rank_genes_group_tracksplot.png')#, save='_{}_rank_genes_group.png'.format(self.name))
                os.chdir("..")

        self.rank_genes_plot_parameter = {'Number of top ranked genes from each group for plot': n_genes}




    ########################### Edit menu stuff ###############################
    def setfigDPI(self):
        wid1 = QInputDialog()
        oldint = self.dpisetting
        temp, okPressed = QInputDialog.getText(wid1, 'DPI setting', 'Set figure DPI:', QLineEdit.Normal,
                                                          str(oldint))
        if (not okPressed) or (not temp.isdigit()):
            self.dpisetting = oldint
            print(temp)
        else:
            self.dpisetting = int(temp)
            if self.whichstage>0:
                sc.settings.set_figure_params(dpi=self.dpisetting, dpi_save=self.dpisetting)
                print("Setting figure DPI to: " + temp)
            else:
                print("Will set figure DPI to: " + temp)

    def clearConsole(self):
        self.plainTextEdit.setPlainText("")

##########################################################################
################ console output stream ##################
    def testlogger(self):
        if self.firstcall:
            self.plainTextEdit.insertPlainText("")
            self.firstcall=False
        else:
            cursor = self.plainTextEdit.textCursor()
            cursor.movePosition(QtGui.QTextCursor.End)
            cursor.insertText("\n")
            #self.plainTextEdit.insertPlainText("\n")
            self.firstcall=True
            self.plainTextEdit.setTextCursor(cursor)
            self.plainTextEdit.ensureCursorVisible()

class XStream(QObject):
    _stdout = None
    _stderr = None

    msgwritten = pyqtSignal(str)

    def flush(self):
        pass

    def fileno(self):
        return -1

    def write(self,msg):
        if not(self.signalsBlocked()):
            self.msgwritten.emit(msg) #np.unicode(msg))

    @staticmethod
    def stdout():
        if (not XStream._stdout):
            XStream._stdout = XStream()
            sys.stdout = XStream._stdout
        return XStream._stdout

    @staticmethod
    def stderr():
        if (not XStream._stderr):
            XStream._stderr = XStream()
            sys.stderr = XStream._stderr
        return XStream._stderr

#########################***************************************************************###########################
#########################**---------------------------*****---------------------------**###########################
#########################***************************************************************###########################



if __name__ == "__main__":
    import sys
    #logging.basicConfig()
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
