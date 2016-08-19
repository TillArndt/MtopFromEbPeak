#!/usr/bin/env python                                                                                                                                                                                       
 
import math, ROOT, json, optparse, os, sys, pprint, array
from ROOT import TCanvas, TGraphErrors, TLine, TF1, TLegend, TPaveText
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
    #input Eb based on input top mass of 166.5 GeV, 169.5 GeV, 171.5, 172.5 GeV, 173.5, 175.5 GeV and 178.5 GeV 
    Eb_input = array('f', [63.90, 65.74, 66.96, 67.57, 68.18, 69.39, 71.20])
    Eb_input_err = array('f', [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    #average Eb as output from the pseudo-experiments
    Eb_output = array('f', [63.17, 64.25, 65.45, 65.95, 66.64, 66.00, 67.99])
    #average uncertainty on Eb as output from the pseudo-experiments
    Eb_output_err = array('f', [0.86, 0.87, 0.83, 0.82, 0.84, 0.86, 0.85])
    #Eb_output_err = array('f', [0.02, 0.02, 0.02])

    graph = TGraphErrors(len(Eb_input), Eb_input, Eb_output, Eb_input_err, Eb_output_err)
    graph.SetMarkerStyle(8)
    graph.SetMarkerSize(0.7)

    fitFunction = TF1("f1", "[0] + [1]*(x-67.57)", 61, 73)
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
    c1.SetLeftMargin(0.15)
    c1.SetBottomMargin(0.15)
    h2 = c1.DrawFrame(61,61,73,73)
    h2.GetXaxis().SetTitle("Predicted energy peak position [GeV]")
    h2.GetXaxis().SetTitleSize(0.045)
    h2.GetXaxis().SetLabelSize(0.06)
    h2.GetXaxis().SetTitleOffset(1.15)
    h2.GetYaxis().SetTitle("Measured energy peak position [GeV]")
    h2.GetYaxis().SetTitleSize(0.045)
    h2.GetYaxis().SetLabelSize(0.06)
    h2.GetYaxis().SetTitleOffset(1.1)
    graph.Draw("p")
    line = TLine(61,61,73,73)
    line.SetLineStyle(2)
    line.SetLineColor(4)
    line.Draw("same")
    fitFunction.SetLineColor(2)
    ROOT.gStyle.SetOptFit(1111)
    ps = graph.GetListOfFunctions().FindObject("stats")
    ps.SetX1NDC(0.55)
    ps.SetY1NDC(0.2)
    ps.SetX2NDC(0.89)
    ps.SetY2NDC(0.35)
    ps.SetBorderSize(0)
    fitFunction.Draw("same")
    leg = TLegend(0.2,0.65,0.5,0.85)
    leg.AddEntry(graph,"Simulation","ep")
    leg.AddEntry(line,"Expected","l")
    leg.AddEntry(fitFunction,"Fit","l")
    leg.SetBorderSize(0)
    leg.Draw("same")
    ##adding text doesn't work yet?
    #cmslabel = TPaveText(0.1,0.2,0.4,0.5)
    #cmslabel.AddText("CMS Simulation")
    #cmslabel.SetFillStyle(0)
    #cmslabel.SetBorderSize(0)
    #cmslabel.SetTextColor(1)
    #cmslabel.SetTextAlign(12)
    #cmslabel.SetTextFont(61)
    #cmslabel.Draw("same")
    c1.Modified()
    c1.Update()
    c1.SaveAs("calibrationCurve.pdf")

    return offset, slope, offset_err, slope_err

def main():

    Ebmeas = 65.73
    offset,slope,offset_err,slope_err = setupCalibrationCurve()
    calibrate(offset,slope,offset_err,slope_err, Ebmeas)

if __name__ == "__main__":
    sys.exit(main())

