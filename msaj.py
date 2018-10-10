#!/usr/bin/env python3

"""
This file is part of MSAJ.

MSAJ is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MSAJ is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MSAJ.  If not, see <https://www.gnu.org/licenses/>.

Copyright 2018 Rodrigo Aluizio
"""

import sys
from PyQt5.QtWidgets import QFileDialog, QAction, QMainWindow, qApp, QTextEdit, QPushButton
from PyQt5.QtGui import QIcon, QTextCursor, QColor
from os.path import join, dirname, basename, splitext
from glob import glob
from Bio import AlignIO
from Bio.Nexus.Nexus import NexusError
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from seq_check import seq_check
from pathlib import Path

if getattr(sys, "frozen", False):
    datadir = dirname(sys.executable)
else:
    datadir = dirname(__file__)


class App(QMainWindow):
 
    def __init__(self):
        super().__init__()

        self.setWindowTitle('Multi Sequence Alignment Joiner')
        self.setGeometry(50, 75, 640, 400)
        self.setFixedSize(640, 400)
        self.setWindowIcon(QIcon(join(datadir, 'icons/dna.svg')))

        self.files = []
        self.msa = []
        self.filename = ''
        self.status = None

        self.open_files = QAction(QIcon(join(datadir, "icons/file.svg")), 'Open Files', self)
        self.open_files.setShortcut('Ctrl+F')
        self.open_files.triggered.connect(self.file_names)
        self.open_files.setEnabled(True)

        self.open_folder = QAction(QIcon(join(datadir, "icons/folder.svg")), 'Open Folder', self)
        self.open_folder.setShortcut('Ctrl+O')
        self.open_folder.triggered.connect(self.folder_name)
        self.open_folder.setEnabled(True)

        self.exit = QAction(QIcon(join(datadir, "icons/exit.svg")), 'Exit', self)
        self.exit.setShortcut('Ctrl+E')
        self.exit.triggered.connect(qApp.quit)
        self.exit.setEnabled(True)

        self.join = QPushButton(' Merge', self)
        self.join.setIcon(QIcon(join(datadir, "icons/merge.svg")))
        self.join.move(275, 355)
        self.join.clicked.connect(self.msaj)
        self.join.setEnabled(False)

        self.toolbar = self.addToolBar('Main')
        self.toolbar.setMovable(False)
        self.toolbar.setStyleSheet('QToolBar{spacing:10px; background:white}')
        self.toolbar.addAction(self.open_files)
        self.toolbar.addAction(self.open_folder)
        self.toolbar.addAction(self.exit)

        self.textbox = QTextEdit(self)
        self.textbox.move(10, 40)
        self.textbox.resize(620, 300)
        self.textbox.setReadOnly(True)

        self.show()

    def file_names(self):
        options = QFileDialog.ReadOnly
        self.files = []
        self.files = QFileDialog.getOpenFileNames(self, 'Select files', str(Path.home()),
                                                  'FASTA Files (*.fas *.fa *.fasta);;'
                                                  'NEXUS Files (*.nex *.nexus *.nxs);;'
                                                  'All Files (*)',
                                                  options=options)[0]
        self.msa = []
        if self.files:
            self.files = sorted(self.files, key=str.lower)
            try:
                for file in self.files:
                    if splitext(basename(file))[1] in ('.fas', '.fa', '.fasta'):
                        self.msa.append(AlignIO.read(file, "fasta", alphabet=IUPACAmbiguousDNA()))
                    if splitext(basename(file))[1] in ('.nex', '.nexus', '.nxs'):
                        self.msa.append(AlignIO.read(file, "nexus", alphabet=IUPACAmbiguousDNA()))
                self.report()
            except (ValueError, NexusError):
                self.textbox.setTextColor(QColor("red"))
                self.textbox.append(f"{basename(file)} is corrupted!")
                self.textbox.setTextColor(QColor("black"))
                self.status = False
        else:
            self.textbox.setTextColor(QColor("red"))
            self.textbox.append("No files were imported!\n")
            self.textbox.setTextColor(QColor("red"))
            self.status = False

    def folder_name(self):
        options = QFileDialog.ReadOnly
        folder = QFileDialog.getExistingDirectory(self, "Select a directory", str(Path.home()), options=options)
        self.msa = []
        self.files = []
        types = ['*.fas', '*.fa', '*.fasta', '*.nex', '*.nexus', '*.nxs']
        for file in types:
            self.files.extend(glob(join(folder, file)))
        if self.files:
            self.files = sorted(self.files, key=str.lower)
            try:
                for file in self.files:
                    if splitext(basename(file))[1] in ('.fas', '.fa', '.fasta'):
                        self.msa.append(AlignIO.read(file, "fasta", alphabet=IUPACAmbiguousDNA()))
                    if splitext(basename(file))[1] in ('.nex', '.nexus', '.nxs'):
                        self.msa.append(AlignIO.read(file, "nexus", alphabet=IUPACAmbiguousDNA()))
                self.report()
            except (ValueError, NexusError):
                self.textbox.setTextColor(QColor("red"))
                self.textbox.append(f"{basename(file)} can not be imported, it may be corrupted or with incorrect"
                                    f" data!")
                self.textbox.setTextColor(QColor("black"))
                self.status = False
        else:
            self.textbox.setTextColor(QColor("red"))
            self.textbox.append("No files were imported!")
            self.textbox.setTextColor(QColor("black"))
            self.status = False

    def report(self):
        diff_ids = seq_check(self.msa)
        if diff_ids == set():
            self.textbox.append("Number of Loci: " + str(len(self.msa)))
            self.textbox.append("Number of Sequences: " + str(len(self.msa[0])))
            self.status = True
            self.join.setEnabled(True)
        else:
            self.textbox.setTextColor(QColor("red"))
            self.textbox.append("Strains: ")
            for value in list(diff_ids)[:-1]:
                self.textbox.insertPlainText(str(value) + ", ")
                self.textbox.moveCursor(QTextCursor.End)
            self.textbox.insertPlainText(str(list(diff_ids)[-1]))
            self.textbox.append("do not have matches in all files!")
            self.textbox.append("Not merging them!")
            self.textbox.setTextColor(QColor("black"))
            self.msa = []
            self.files = []
            self.status = False

    def save_file(self):
        options = QFileDialog.Options()
        self.filename = QFileDialog.getSaveFileName(self, "Save File", str(Path.home()),
                                                    "NEXUS Files (*.nex *.nexus *.nxs)", options=options)[0]
        if basename(splitext(self.filename)[1]) == '' and basename(splitext(self.filename)[0]) != '':
            self.filename += '.nex'
        return self.filename

    def msaj(self):
        self.msa[0].sort()
        combined = self.msa[0]
        for i in range(1, len(self.msa)):
            self.msa[i].sort()
            combined = combined + self.msa[i]

        # Save multi locus multiple sequence alignment
        try:
            with open(self.save_file(), "w") as f:
                AlignIO.write(combined, f, "nexus")
                f.write('\nbegin sets;\n')
                count = 1
                for i in range(len(self.files)):
                    f.write(f'\tcharset {splitext(basename(self.files[i]))[0]} = {str(count)}-'
                            f'{str(self.msa[i].get_alignment_length() + count - 1)}\n')
                    count = count + self.msa[i].get_alignment_length()
                f.write('end;\n')

            # Print some simple but useful information
            self.textbox.append("\nMulti Locus Multiple Sequence Alignment created")
            self.textbox.append(f"at {self.filename}")
            self.textbox.append("\n----- Useful Information -----")
            count = 1
            for i in range(len(self.files)):
                self.textbox.append("Locus " + splitext(basename(self.files[i]))[0] + ": " + str(count) +
                                    " - " + str(self.msa[i].get_alignment_length() + count - 1))
                count = count + self.msa[i].get_alignment_length()
            self.textbox.append("-----------------------------------------\n")

            self.msa = []
            self.files = []
            self.join.setEnabled(False)
        except (TypeError, IndexError, FileNotFoundError):
            self.msa = []
            self.files = []
            self.join.setEnabled(False)
