# HVM3 Calibration Software

You can calibrate an HVM3 image with the following syntax:

> python hvm3_calibrate.py file_list.txt

... where file_list.txt is an ASCII text file in the following format:

<dark_data_file_1> <science_data_file_1>

<dark_data_file_2> <science_data_file_2>

<dark_data_file_3> <science_data_file_3>

...et cetera.

Directory contents:
- scripts/  contains scripts used to generate the calibration files.
- config/  contains calibration configuration information.
- data/ contains the calibration data files
- utils/ contains executable utilities and library functions used in calibration.

----------------
Copyright (c) 2023-24 California Institute of Technology (“Caltech”). U.S. Government sponsorship acknowledged.
All rights reserved.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
