#!/usr/bin/env python                                                                                                                                                                                       
 
import math, ROOT, json, optparse, os, sys, pprint, array
from ROOT import TCanvas, TGraphErrors, TLine, TF1
from array import array
from math import *

def calibrate(offset=None,slope=None,offset_err=None,slope_err=None,meas=None):
    
    #first get calibrated Eb from the calibration curve and the measured value
    calibrated = (meas - offset)/slope + 67.57
    m_W = 80.385
    m_b = 4.18
    m_t = calibrated + sqrt(m_W*m_W - m_b*m_b + calibrated*calibrated)
    return m_t

def setupCalibrationCurve():
    #input Eb based on input top mass of 169.5 GeV, 171.5 GeV, 172.5 GeV, 173.5 GeV, 175.5 GeV 
    Eb_input = array('f', [65.74, 67.57, 69.39])
    Eb_input_err = array('f', [0.0, 0.0, 0.0])
    #average Eb as output from the pseudo-experiments
    Eb_output = array('f', [64.74, 66.81, 67.85])
    #average uncertainty on Eb as output from the pseudo-experiments
    Eb_output_err = array('f', [0.69, 0.76, 0.74])
    #Eb_output_err = array('f', [0.02, 0.02, 0.02])

    graph = TGraphErrors(len(Eb_input), Eb_input, Eb_output, Eb_input_err, Eb_output_err)
    graph.SetMarkerStyle(8)
    graph.SetMarkerSize(0.7)

    fitFunction = TF1("f1", "[0] + [1]*(x-67.57)", 64, 70)
#    fitFunction = TF1("f1", "[0] + [1]*(x)", 64, 70)
    fitFunction.SetLineColor(3)
    fitRes = graph.Fit("f1","S")
    offset = fitRes.Parameter(0)
    offset_err = fitRes.ParError(0)
    slope = fitRes.Parameter(1)
    slope_err = fitRes.ParError(1)
    chi = fitRes.Chi2()
    cor = fitRes.Correlation(0,1)

    print "offset = (%3.2f #pm %3.2f)" % (offset, offset_err)
    print "slope = (%3.2f #pm %3.2f)" % (slope, slope_err)
    print "chi2 = %3.2f" % chi
    print "correlation = %3.2f" % cor

    c1 = TCanvas("c1", "c1",0,0,600,500);
    c1.Range(0,0,1,1)
    h2 = c1.DrawFrame(64,64,70,70)
    h2.GetXaxis().SetTitle("Eb_{input}")
    h2.GetXaxis().SetTitleSize(0.06)
    h2.GetXaxis().SetLabelSize(0.06)
    h2.GetXaxis().SetTitleOffset(0.45)
    h2.GetYaxis().SetTitle("Eb_{output}")
    h2.GetYaxis().SetTitleSize(0.06)
    h2.GetYaxis().SetLabelSize(0.06)
    h2.GetYaxis().SetTitleOffset(0.45)
    graph.Draw("p")
    line = TLine(64,64,70,70)
    line.Draw("same")
    fitFunction.Draw("same")
    c1.Modified()
    c1.Update()
    c1.SaveAs("calibrationCurve.pdf")

    return offset, slope, offset_err, slope_err

def main():

    offset,slope,offset_err,slope_err = setupCalibrationCurve()
    #test
    calibrate(offset,slope,offset_err,slope_err, 67.02)

if __name__ == "__main__":
    sys.exit(main())

