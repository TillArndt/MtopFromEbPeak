#!/usr/bin/env python                                                                                                                                                                                       
 
import math, ROOT, json, optparse, os, sys, pprint, array
from ROOT import TCanvas, TGraphErrors, TLine, TF1
from array import array
from math import *

def main():

    m_top = 172.5
    m_w = 80.385
    m_b = 4.18
    Eb = math.sqrt( (m_top*m_top + m_w*m_w - m_b*m_b) * (m_top*m_top + m_w*m_w - m_b*m_b) / ((2 * m_top) * (2 * m_top)) - m_w*m_w + m_b*m_b);
    print "predicted Eb for mass %3.2f GeV = %3.2f GeV" % (m_top,Eb)

if __name__ == "__main__":
    sys.exit(main())
