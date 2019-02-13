#ifndef PLOTMAKER_H_
#define PLOTMAKER_H_

#include "TString.h"
#include "TList.h"
#include "THStack.h"

#include "RooWorkspace.h"
#include "RooCategory.h"

#include "workspace_helper.h"

class LabelInfo: public TObject {
    public:
        TString name;
        TString label;
        Color_t clr;

        LabelInfo(TString n, TString l, Color_t c): name(n), label(l), clr(c) {}
        ~LabelInfo() {}

        ClassDef(LabelInfo,1);
};

class PartitionInfo: public TObject {
    public:
        TString name;
        Int_t bin_lo;
        Int_t bin_hi;
        Int_t nbins;

        PartitionInfo(TString n, Int_t b1, Int_t b2, Int_t bins): name(n), bin_lo(b1), bin_hi(b2), nbins(bins) {}
        ~PartitionInfo() {}

        ClassDef(PartitionInfo,1);
};

class PlotMaker {
    private:
        typedef std::vector<TString> vTStr;
        typedef std::vector<RooAddition*> vRooAdd;
        typedef std::pair<TString,double> pTStrDbl;
        typedef PartitionInfo PInfo;

    public:
        PlotMaker(RooWorkspace* ws,TString out_dir,TString ftype);
        ~PlotMaker();

        void setFitResult(RooFitResult* full, RooFitResult* stat);
        void printObsData(std::vector<pTStrDbl> data, vTStr cats);
        void printSummedYields(vRooAdd summed_yields,vTStr cats);
        void makeYieldsPlot(
            std::vector<pTStrDbl> data,
            vRooAdd yields,
            vTStr procs,
            TString fname,
            TString plot_title=""
        );
        void makeStackPlot(
            TString fname,
            TH1D* h_data,
            TH1* ratio,
            TH1* r_1sig,
            TH1* r_2sig,
            TH1* h_err,
            TH1* h_stats,
            std::vector<TH1*> hists,
            std::vector<PInfo> partitions={}
        );
        void makeSummedPlot(
            TString fname,
            TH1D* h_data,
            TH1* ratio,
            TH1* r_1sig,
            TH1* r_2sig,
            std::vector<TH1D*> hists
        );
        void makeCorrelationPlot(TString fname,RooFitResult* fr);
        void fillDataHist(
            TH1D* h_data,
            std::vector<pTStrDbl> data,
            vTStr cats
        );
        void fillRatioHists(
            TH1D* h_data,
            TH1* ratio,
            TH1* r_1sig,
            TH1* r_2sig,
            TH1* h_err,
            vRooAdd yields,
            vTStr cats
        );
        void fillErrorHist(TH1* h_err, TH1* h_stats, vRooAdd yields, vTStr cats);
        void fillProcessYieldHists(
            TH1D* h_data,
            std::vector<TH1*> *hists,
            vRooAdd yields,
            vTStr cats,
            vTStr procs
        );

        int findCatIndex(vTStr cats, TString search);
        void setupCanvas(TCanvas* canv, bool inc_ratio=true);
        void setupHists(TH1D* data, TH1* ratio, TH1* ratio_1sig, TH1* ratio_2sig);
        void setupLegends(TLegend *leg, TLegend *leg_ratio);

        TString output_dir;
        TString output_type;
        RooWorkspace* ws;
        std::map<TString,LabelInfo*> proc_map;

        bool nodata;
        bool noratio;
        bool isfakedata;
        bool splitunc;

        WSHelper ws_helper;

        RooFitResult* fr_full;
        RooFitResult* fr_stat;
};

PlotMaker::PlotMaker(RooWorkspace* ws,TString out_dir,TString ftype) {
    this->ws = ws;
    this->output_dir  = out_dir;
    this->output_type = ftype;
    this->ws_helper   = WSHelper();

    this->proc_map["ttH"] = new LabelInfo("ttH","t#bar{t}H",kRed+1);
    this->proc_map["ttll"] = new LabelInfo("ttll","ttll",kGreen+1);
    this->proc_map["ttlnu"] = new LabelInfo("ttlnu","ttl#nu",kBlue);
    this->proc_map["tllq"] = new LabelInfo("tllq","tllq",kPink+1);
    this->proc_map["charge_flips"] = new LabelInfo("charge_flips","Flips",kAzure-9);
    this->proc_map["fakes"] = new LabelInfo("fakes","Fakes",kYellow-7);
    this->proc_map["WZ"] = new LabelInfo("WZ","WZ",kMagenta);
    this->proc_map["ZZ"] = new LabelInfo("ZZ","ZZ",kViolet+1);
    this->proc_map["WW"] = new LabelInfo("WW","WW",kViolet+2);
    this->proc_map["WWW"] = new LabelInfo("WWW","WWW",kSpring+1);
    this->proc_map["WWZ"] = new LabelInfo("WWZ","WWZ",kSpring+2);
    this->proc_map["WZZ"] = new LabelInfo("WZZ","WZZ",kSpring+3);
    this->proc_map["ZZZ"] = new LabelInfo("ZZZ","ZZZ",kSpring+4);
    this->proc_map["ttGJets"] = new LabelInfo("convs","Convs",kYellow+3);
    this->proc_map["singlet_tWchan"] = new LabelInfo("tW","tW",kMagenta+1);
    this->proc_map["singletbar_tWchan"] = new LabelInfo("tbarW","#bar{t}W",kMagenta+2);
    this->proc_map["diboson"] = new LabelInfo("diboson","Diboson",kMagenta);
    this->proc_map["triboson"] = new LabelInfo("triboson","Triboson",kSpring+1);
    this->proc_map["MC"] = new LabelInfo("MC","MC",kGray);

    this->nodata = true;
    this->noratio = false;
    this->isfakedata = true;
    this->splitunc = true;      // Split uncertainty into syst and stat err
}

PlotMaker::~PlotMaker() {}

void PlotMaker::setFitResult(RooFitResult* full, RooFitResult* stat) {
    this->fr_full = full;
    this->fr_stat = stat;
}

void PlotMaker::printObsData(std::vector<pTStrDbl> data, vTStr cats) {
    // Print out out the observed yields in each category
    TString indent = "    ";

    std::cout << "Obs Data:" << std::endl;

    int dMax = 0;
    for (auto c: cats) {
        dMax = (dMax > c.Length()) ? dMax : c.Length();
    }

    for (auto c: cats) {
        for (auto d: data) {
            if (c == d.first) {
                std::cout << indent << std::left << std::setw(dMax) << d.first << ": " << d.second << std::endl;
                break;
            }
        }
    }
}

void PlotMaker::printSummedYields(vRooAdd summed_yields,vTStr cats) {
    // Prints the yields for categories which have been summed over certain processes
    TString indent = "    ";

    std::cout << "Sum Yields:" << std::endl;
    if (summed_yields.size() != cats.size()) {
        std::cout << indent << "Size mis-match between summed_yields and categories!" << std::endl;
        return;
    }
    int dMax = 0;
    for (auto c: cats) {
        dMax = (dMax > c.Length()) ? dMax : c.Length();
    }
    for (uint i = 0; i < summed_yields.size(); i++) {
        RooAddition* summed_cat = summed_yields.at(i);
        TString c = cats.at(i);

        double val = summed_cat->getVal();
        double full_err = summed_cat->getPropagatedError(*(this->fr_full));
        double stat_err = summed_cat->getPropagatedError(*(this->fr_stat));
        double syst_err = sqrt(pow(full_err,2) - pow(stat_err,2));

        if (this->splitunc) {
            std::cout << indent << std::left << std::setw(dMax) << c << ": "
                      << std::setw(8) << summed_cat->getVal() << " +/- "
                      << std::setw(9) << full_err << " | "
                      << std::setw(9) << stat_err << " (stat) " << std::setw(9) << syst_err << " (syst)" << std::endl;
        } else {
            std::cout << indent << std::left << std::setw(dMax) << c << ": "
                      << std::setw(8) << summed_cat->getVal() << " +/- "
                      << std::setw(9) << full_err << std::endl;
        }
    }
}

void PlotMaker::makeYieldsPlot(
    std::vector<pTStrDbl> data,
    vRooAdd yields,
    vTStr procs,
    TString fname,
    TString plot_title=""
) {
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1,0);

    TString indent = "    ";
    TString line_break = "*************************************************************";

    //std::cout << "All Yields:" << std::endl;
    for (RooAddition* f: yields) {// Remove the 'n_exp_bin' from the naming if needed
        TString old_name = f->GetName();
        TString new_name = old_name.ReplaceAll("n_exp_bin","");
        f->SetName(new_name);
        //std::cout << indent << f->GetName() << ": " << f->getVal() << std::endl;
    }

    vTStr cats = this->ws_helper.stripProcessName(yields);   // This determines the bin ordering!
    int ncats = cats.size();

    this->printObsData(data,cats);

    vRooAdd summed_yields;
    for (auto c: cats) {// Sum over all MC processes
        RooAddition* summed_cat = this->ws_helper.sumProcesses(yields,c,procs);
        summed_yields.push_back(summed_cat);
    }

    this->printSummedYields(summed_yields,cats);

    //return;

    std::vector<TString> bkg_procs {"charge_flips","fakes","ttGJets","WZ","ZZ","WW","WWW","WWZ","WZZ","ZZZ",
                                    "singlet_tWchan","singletbar_tWchan"};
    vRooAdd summed_bkg;
    for (auto c: cats) {
        RooAddition* summed_cat = this->ws_helper.sumProcesses(yields,c,bkg_procs);
        summed_bkg.push_back(summed_cat);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////

    std::vector<TH1*> proc_yields;
    TH1D* h_data = new TH1D("data_yield","data;category",ncats,0.0,1.0);//ncats

    TH1* h_err   = (TH1*)h_data->Clone("err_hist");
    TH1* ratio   = (TH1*)h_data->Clone("ratio");
    TH1* r_1sig  = (TH1*)ratio->Clone("ratio_1sig");
    TH1* r_2sig  = (TH1*)ratio->Clone("ratio_2sig");

    h_err->Reset();
    h_err->SetFillColor(1);
    h_err->SetFillStyle(3244);//Criss-crossed hashing
    h_err->SetLineColor(1);
    h_err->SetMarkerColor(0);
    h_err->SetMarkerStyle(1);
    h_err->SetMarkerSize(0.);

    TH1* h_stats = (TH1*)h_err->Clone("stat_err");

    if (this->splitunc) {
        h_err->SetFillStyle(3445);
        h_stats->SetFillStyle(3244);
    } else {
        /*  Styles:
                3345 - Single diagonal hash
                1001 - Solid
        */
        h_err->SetFillStyle(3345);
        //h_err->SetFillStyle(1001); h_err->SetFillColor(kGray); h_err->SetFillColorAlpha(kGray,0.5);
    }

    //if (this->nodata) {
    //    // TODO: Switch this to the full error histogram
    //    h_err->SetFillStyle(3244);
    //    h_stats->SetFillStyle(3244);
    //}

    this->setupHists(h_data,ratio,r_1sig,r_2sig);

    if (plot_title.Length()) {// Set customized title
        h_data->SetTitle(plot_title+";;Events");
    }

    this->fillDataHist(h_data,data,cats);
    this->fillErrorHist(h_err,h_stats,summed_yields,cats);
    this->fillRatioHists(h_data,ratio,r_1sig,r_2sig,h_err,summed_yields,cats);
    this->fillProcessYieldHists(h_data,&proc_yields,yields,cats,procs);

    TH1D* summed_hist = (TH1D*)h_data->Clone("h_summed_ex");
    summed_hist->Reset();

    //for (auto h: proc_yields) {
    //    summed_hist->Add(h,1.0);
    //}

    /*
    std::cout << "Proc Yield Val:" << std::endl;
    std::cout << indent << std::setw(13) << "";
    for (uint i = 0; i < procs.size(); i++) {
        TString p = procs.at(i);
        if (p.Length() >= 8) {
            p = p(0,8);
        }
        if (i) {
            std::cout << " + " << std::setw(8) << p;
        } else {
            std::cout << std::setw(8) << p;
        }
    }
    std::cout << std::endl;
    for (auto c: cats) {
        int cat_idx = this->findCatIndex(cats,c); // Could also just use a standard 'for' loop instead
        int bin_idx = cat_idx + 1;
        std::cout << indent << std::setw(13) << (c+":");
        for (uint i = 0; i < proc_yields.size(); i++) {
            double bin_err = proc_yields.at(i)->GetBinContent(bin_idx);
            if (i) {
                std::cout << " + " << std::setw(8) << bin_err;
            } else {
                std::cout << std::setw(8) << bin_err;
            }
        }
        std::cout << std::endl;
    }

    std::cout << "Proc Yield Err:" << std::endl;
    for (auto c: cats) {
        int cat_idx = this->findCatIndex(cats,c); // Could also just use a standard 'for' loop instead
        int bin_idx = cat_idx + 1;
        std::cout << indent << std::setw(13) << (c+":");
        for (uint i = 0; i < proc_yields.size(); i++) {
            double bin_err = proc_yields.at(i)->GetBinError(bin_idx);
            if (i) {
                std::cout << " + " << std::setw(8) << bin_err;
            } else {
                std::cout << std::setw(8) << bin_err;
            }
        }
        std::cout << std::endl;
    }

    std::cout << "Summed Plot Err:" << std::endl;
    for (auto c: cats) {
        int cat_idx = this->findCatIndex(cats,c); // Could also just use a standard 'for' loop instead
        int bin_idx = cat_idx + 1;
        double bin_err = summed_hist->GetBinError(bin_idx);
        std::cout << indent << std::setw(13) << (c+":") << " " << bin_err << std::endl;
    }
    */

    std::vector<PInfo> partitions;

    //PInfo p_ee   = PInfo("e^{+}e^{+}",0,5,34);
    //PInfo p_emu  = PInfo("e^{+}#mu^{+}",5,12,34);
    //PInfo p_mumu = PInfo("#mu^{+}#mu^{+}",12,17,34);
    //PInfo m_ee   = PInfo("e^{-}e^{-}",17,24,34);
    //PInfo m_emu  = PInfo("e^{-}#mu^{-}",24,31,34);
    //PInfo m_mumu = PInfo("#mu^{-}#mu^{-}",31,34,34);

    //partitions.push_back(p_ee);
    //partitions.push_back(p_emu);
    //partitions.push_back(p_mumu);
    //partitions.push_back(m_ee);
    //partitions.push_back(m_emu);
    //partitions.push_back(m_mumu);

    TH1* empty_ratio = 0;
    this->makeStackPlot(fname + "_stack",h_data,ratio,r_1sig,r_2sig,h_err,h_stats,proc_yields,partitions);
    //this->makeSummedPlot(fname + "_sum",h_data,ratio,r_1sig,r_2sig,proc_yields);

    //delete h_data;
    //delete ratio;
    //delete r_1sig;
    //delete r_2sig;
}

void PlotMaker::makeStackPlot(
    TString fname,
    TH1D* h_data,
    TH1* ratio,TH1* r_1sig,TH1* r_2sig, TH1* h_err, TH1* h_stats,
    std::vector<TH1*> hists,
    std::vector<PInfo> partitions={}
) {
    bool include_ratio = true;
    if (!ratio || this->noratio || this->nodata) {
        include_ratio = false;
    }

    TCanvas* canv = new TCanvas("canv","canv",960,700);//(960,943) or (600,700)
    this->setupCanvas(canv,include_ratio);

    TLegend *d_leg = new TLegend(0.14,0.75,0.94,0.89); //(x1,y1,x2,y2) (0.20,0.80,0.94,0.89)
    TLegend *r_leg = new TLegend(0.14,0.77,0.30,0.97); //(x1,y1,x2,y2) (0.14,0.77,0.48,0.97)
    this->setupLegends(d_leg,r_leg);

    //if (!(this->nodata)) d_leg->AddEntry(h_data,"data","lp");//lpe
    //d_leg->AddEntry(h_err,"Stat+Syst","f");
    //d_leg->AddEntry(h_stats,"Stat","f");
    r_leg->AddEntry(r_1sig,"Stat+Syst","f");
    r_leg->AddEntry(r_2sig,"Stat","f");

    if (this->nodata) {
        d_leg->AddEntry(h_err,"Total unc.","f");
    } else {
        if (this->isfakedata) {
            d_leg->AddEntry(h_data,"fake data","lp");
        } else {
            d_leg->AddEntry(h_data,"data","lp");
        }
        if (this->splitunc) {
            d_leg->AddEntry(h_err,"Stat+Syst","f");
            d_leg->AddEntry(h_stats,"Stat","f");
        } else {
            d_leg->AddEntry(h_err,"Total unc.","f");
        }
    }

    THStack *hs = new THStack("hs_category_yield","");
    for (auto h: hists) {
        TString p = h->GetName();
        d_leg->AddEntry(h,this->proc_map[p]->label,"f");
        hs->Add(h);
    }

    int max_data_bin = h_data->GetMaximumBin();
    double max_data  = h_data->GetBinContent(max_data_bin) + h_data->GetBinError(max_data_bin);
    double max_stack = hs->GetMaximum();
    double hmax = this->nodata ? max_stack : TMath::Max(max_data,max_stack);
    //h_data->SetMaximum(hmax*1.5);
    //h_data->SetMinimum(0.1);

    hs->SetMaximum(hmax*1.5);
    hs->SetMinimum(0.0);//0.1
    //hs->SetTitle(";;Events");
    TString title = h_data->GetTitle();
    title += ";;Events";
    hs->SetTitle(title);

    canv->cd(1);

    hs->Draw("hist");
    hs->GetYaxis()->SetTitleSize(0.05);
    hs->GetYaxis()->SetTitleOffset(1.1);
    hs->GetYaxis()->SetLabelSize(0.05);

    // I have no idea why this is the only place where it is able re-align the tick marks
    hs->GetXaxis()->SetNdivisions(-1,kFALSE);

    h_err->Draw("e2,same");
    if (this->splitunc) {
        h_stats->Draw("e2,same");
    }

    if (!(this->nodata)) h_data->Draw("e,l,same");

    /*
    //h_data->Draw("e,l");
    //hs->Draw("hist,same");
    hs->Draw("hist");
    h_data->Draw("e,l,same");
    */

    double line_height = h_data->GetMaximum()*0.8;
    for (auto p: partitions) {
        Int_t p_width = p.bin_hi - p.bin_lo;
        TLine* l = new TLine(p.bin_hi + 1,0.0,p.bin_hi + 1,line_height);
        l->SetLineColor(kRed);
        l->SetLineWidth(2);
        l->SetLineStyle(kSolid);//kSolid = 1, kDashed, kDotted, kDashDotted
        if (p.bin_hi != p.nbins) {
            l->Draw();
        }

        double left_margin = .11;
        double right_margin = .05;

        double width = (1.0 - left_margin - right_margin) / p.nbins;
        double mid_pt = left_margin + (p.bin_lo + p.bin_hi)*width/2.0;
        //double mid_pt = (p.bin_lo + p.bin_hi)/(2.0);


        std::cout << "Point: " << mid_pt << std::endl;
        TLatex* latex_label = new TLatex(mid_pt,0.02,p.name);
        latex_label->SetNDC();
        latex_label->SetTextFont(42);
        latex_label->SetTextSize(0.035);//0.05
        latex_label->SetTextAlign(21);
        latex_label->Draw();
    }

    d_leg->Draw();

    if (include_ratio) {
        canv->cd(2);
        r_1sig->Draw("e2");
        r_2sig->Draw("e2,same");
        ratio->Draw("e1,same");

        canv->GetPad(1)->RedrawAxis();
        canv->GetPad(2)->RedrawAxis();
        
        //r_leg->Draw();
    } else {
        canv->GetPad(1)->RedrawAxis();
    }

    TString save_name = this->output_dir + "/" + fname + "." + this->output_type;
    canv->Print(save_name,this->output_type);

    //delete hs;
    //delete canv;
    //delete d_leg;
    //delete r_leg;
}

void PlotMaker::makeSummedPlot(
    TString fname,
    TH1D* h_data,
    TH1* ratio, TH1* r_1sig, TH1* r_2sig,
    std::vector<TH1D*> hists
) {
    TCanvas* canv = new TCanvas("canv","canv",960,700);//(960,943) or (600,700)
    this->setupCanvas(canv);

    TLegend *d_leg = new TLegend(0.14,0.75,0.94,0.89); //(x1,y1,x2,y2) (0.2,0.8,0.94,0.89)
    TLegend *r_leg = new TLegend(0.14,0.77,0.48,0.97);
    this->setupLegends(d_leg,r_leg);

    d_leg->AddEntry(h_data,"data","lpe");
    r_leg->AddEntry(r_1sig,"Stat Uncertainty","f");

    TString p = "MC";
    TH1D* hist = (TH1D*)h_data->Clone("h_category_summed_yield");
    hist->Reset();

    //hist->SetFillColor(kWhite);
    hist->SetFillColor(this->proc_map[p]->clr);
    hist->SetLineColor(this->proc_map[p]->clr);

    for (TH1D* h: hists) {
        hist->Add(h,1.0);
    }

    d_leg->AddEntry(hist,this->proc_map[p]->label,"f");

    canv->cd(1);

    h_data->Draw("e,l");
    hist->Draw("hist,e,same");
    h_data->Draw("e,l,same");

    d_leg->Draw();

    canv->cd(2);
    r_1sig->Draw("e2");
    ratio->Draw("e1,same");

    canv->GetPad(1)->RedrawAxis();
    canv->GetPad(2)->RedrawAxis();

    r_leg->Draw();

    TString save_name = this->output_dir + "/" + fname + "." + this->output_type;
    canv->Print(save_name,this->output_type);

    delete hist;
    delete canv;
    delete d_leg;
    delete r_leg;
}

void PlotMaker::makeCorrelationPlot(TString fname,RooFitResult* fr) {
    TCanvas* c = new TCanvas("canv","canv",960,700);

    TH2D* h2 = (TH2D*)fr->correlationHist();

    h2->Draw("colz");

    TString save_name = this->output_dir + "/" + fname + "." + this->output_type;
    c->Print(save_name,this->output_type);

    delete h2;
    delete c;
}

// Fill in the data histogram
void PlotMaker::fillDataHist(TH1D* h_data, std::vector<pTStrDbl> data, vTStr cats) {
    for (auto cat: cats) {
        int cat_idx = this->findCatIndex(cats,cat); // Could also just use a standard 'for' loop instead
        int bin_idx = cat_idx + 1;
        for (auto d: data) {
            if (cat == d.first) {
                h_data->SetBinContent(bin_idx,d.second);
                h_data->SetBinError(bin_idx,0.001);//h_data->SetBinError(bin_idx,sqrt(d.second));
                break;
            }
        }
    }
}

// Fill in the ratio histograms (also sets the bin labels)
void PlotMaker::fillRatioHists(TH1D* h_data, TH1* ratio, TH1* r_1sig, TH1* r_2sig, TH1* h_err, vRooAdd yields, vTStr cats) {
    double r_min = 0.0; // Maybe make these data members of PlotMaker
    double r_max = 2.3; // Maybe make these data members of PlotMaker

    r_1sig->SetMinimum(r_min);
    r_1sig->SetMaximum(r_max);

    for (auto cat: cats) {
        int cat_idx = this->findCatIndex(cats,cat);
        int bin_idx = cat_idx + 1;  // Histogram bins are offset by 1, since idx 0 is underflow bin

        double data_central = h_data->GetBinContent(bin_idx);
        double mc_central = yields.at(cat_idx)->getVal();

        double mc_err  = yields.at(cat_idx)->getPropagatedError(*(this->fr_full));
        double mc_stat = yields.at(cat_idx)->getPropagatedError(*(this->fr_stat));
        double mc_syst = sqrt(pow(mc_err,2) - pow(mc_stat,2));

        //double mc_up   = mc_central + mc_err/2.0;
        //double mc_down = mc_central - mc_err/2.0;
        //double ratio_up   = mc_up / mc_central;
        //double ratio_down = mc_down / mc_central;

        double ratio_val = data_central / mc_central;

        //double env_err = sqrt(data_central) / mc_central;
        //double env_err = ratio_val*sqrt(pow(mc_err/mc_central,2));//ratio_val*mc_err/mc_central;
        //double env_err = data_central*mc_err/pow(mc_central,2);
        //double env_err = ratio_up - ratio_down;
        

        //double env_err  = mc_err / mc_central;
        //double stat_err = mc_stat / mc_central;

        double env_err  = data_central*mc_err/(mc_central*mc_central);
        double stat_err = data_central*mc_stat/(mc_central*mc_central);

        ratio->SetBinContent(bin_idx,ratio_val);
        ratio->SetBinError(bin_idx,0.001);

        r_1sig->SetBinContent(bin_idx,ratio_val);
        r_1sig->SetBinError(bin_idx,env_err);

        r_2sig->SetBinContent(bin_idx,ratio_val);
        r_2sig->SetBinError(bin_idx,stat_err);

        if (ratio_val > r_max && (ratio_val - env_err) < r_max) {
            double minner = ratio_val - env_err;
            ratio->SetBinContent(bin_idx,r_max-0.0001);
            ratio->SetBinError(bin_idx,r_max-0.0001-minner);
        }

        TString axis_label = cat.ReplaceAll("C_2lss_","");
        axis_label = axis_label.ReplaceAll("C_3l_","");
        axis_label = axis_label.ReplaceAll("C_4l_","");

        //axis_label = axis_label.ReplaceAll("mix_sfz_","");
        //axis_label = axis_label.ReplaceAll("mix_","");

        //axis_label = axis_label.ReplaceAll("p_ee_","");
        //axis_label = axis_label.ReplaceAll("p_emu_","");
        //axis_label = axis_label.ReplaceAll("p_mumu_","");
        //axis_label = axis_label.ReplaceAll("m_ee_","");
        //axis_label = axis_label.ReplaceAll("m_emu_","");
        //axis_label = axis_label.ReplaceAll("m_mumu_","");

        r_1sig->GetXaxis()->SetBinLabel(bin_idx,axis_label);   // Label the x-axis with categories
        h_data->GetXaxis()->SetBinLabel(bin_idx,axis_label);
    }   
}

void PlotMaker::fillErrorHist(TH1* h_err, TH1* h_stats, vRooAdd yields, vTStr cats) {
    for (auto cat: cats) {
        int cat_idx = this->findCatIndex(cats,cat);
        int bin_idx = cat_idx + 1;  // Histogram bins are offset by 1, since idx 0 is underflow bin

        double mc_central = yields.at(cat_idx)->getVal();
        double full_err = yields.at(cat_idx)->getPropagatedError(*(this->fr_full));
        double stat_err = yields.at(cat_idx)->getPropagatedError(*(this->fr_stat));
        double syst_err = sqrt(pow(full_err,2) - pow(stat_err,2));

        h_err->SetBinContent(bin_idx,mc_central);
        h_err->SetBinError(bin_idx,full_err);

        h_stats->SetBinContent(bin_idx,mc_central);
        h_stats->SetBinError(bin_idx,stat_err);
    }
}

// Create and fill yield histogram for each specified process
void PlotMaker::fillProcessYieldHists(TH1D* clone_hist, std::vector<TH1*> *hists, vRooAdd yields, vTStr cats, vTStr procs) {
    TString col_sep = " | ";
    for (TString p: procs) {
        //std::cout << "Filling Yield: " << p << std::endl;
        //TH1D* hist = (TH1D*)clone_hist->Clone(p);
        TH1* hist = (TH1*)clone_hist->Clone(p);
        hist->Reset();
        hist->SetFillColor(this->proc_map[p]->clr);
        hist->SetLineColor(this->proc_map[p]->clr);
        hist->SetLineWidth(0);
        vRooAdd exp_proc_cats = this->ws_helper.filter(yields,"proc_"+p+"$");
        //std::cout << std::left << std::setw(4) << "BIdx" << col_sep
        //          << std::left << std::setw(22) << "Bin Name" << col_sep
        //          << std::left << std::setw(10) << "Center" << col_sep
        //          << std::left << std::setw(10) << "Low Edge" << col_sep
        //          << std::left << std::setw(10) << "High Edge" << col_sep
        //          << std::left << std::setw(10) << "Width" << col_sep
        //          << std::left << std::setw(10) << "Value" << col_sep
        //          << std::endl;
        for (RooAddition* exp_yield: exp_proc_cats) {
            TString name = exp_yield->GetName();
            //std::cout << "  " << std::setw(18+p.Length()) << (name+":") << " " << exp_yield->getVal() << std::endl;
            TString cat = this->ws_helper.stripProcessName(exp_yield);
            int cat_idx = this->findCatIndex(cats,cat);
            int bin_idx = cat_idx + 1;
            double bin_val = exp_yield->getVal();
            double bin_err = sqrt(bin_val);

            double bin_center = hist->GetBinCenter(bin_idx);
            double bin_width = hist->GetBinWidth(bin_idx);

            double low_edge = hist->GetXaxis()->GetBinLowEdge(bin_idx);
            double high_edge = hist->GetXaxis()->GetBinUpEdge(bin_idx);

            hist->Fill(bin_center,bin_val);

            //hist->Fill(cat_idx,bin_val);

            //hist->SetBinContent(bin_idx,bin_val);
            //hist->SetBinError(bin_idx,bin_err);

            //std::cout << std::left << std::setw(4) << bin_idx << col_sep
            //          << std::left << std::setw(22) << cat << col_sep
            //          << std::left << std::setw(10) << bin_center << col_sep
            //          << std::left << std::setw(10) << low_edge << col_sep
            //          << std::left << std::setw(10) << high_edge << col_sep
            //          << std::left << std::setw(10) << bin_width << col_sep
            //          << std::left << std::setw(10) << bin_val << col_sep
            //          << std::endl;

        }
        hists->push_back(hist);
    }
}

// Attempt to find the index of the specified category
int PlotMaker::findCatIndex(vTStr cats, TString search) {
    int idx = 0;
    for (auto s: cats) {
        if (s == search) return idx;
        idx++;
    }
    return -1;
}

// Basic settings for a canvas to display a histogram with a ratio plot
void PlotMaker::setupCanvas(TCanvas* canv,bool inc_ratio=true) {
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    Float_t small = 1.e-5;
    const float padding = 1e-5;
    const float ydiv = 0.3;
    if (inc_ratio) {    
        canv->Divide(1,2,small,small);
        canv->GetPad(1)->SetPad(padding,ydiv+padding,1-padding,1-padding);
        canv->GetPad(1)->SetLeftMargin(.11);
        canv->GetPad(1)->SetRightMargin(.05);
        canv->GetPad(1)->SetBottomMargin(.3);
        canv->GetPad(1)->Modified();

        canv->GetPad(2)->SetPad(padding,padding,1-padding,ydiv-padding);
        canv->GetPad(2)->SetLeftMargin(.11);
        canv->GetPad(2)->SetRightMargin(.05);
        canv->GetPad(2)->SetBottomMargin(.3);
        canv->GetPad(2)->Modified();

        canv->cd(1);
        gPad->SetBottomMargin(small);
        gPad->Modified();

        canv->cd(2);
        gPad->SetTopMargin(small);
        gPad->SetTickx();
        gPad->Modified();
    } else {
        canv->Divide(1,1,small,small);
        //canv->GetPad(1)->SetPad(xlow,ylow,xup,yup);
        canv->GetPad(1)->SetPad(padding,padding,1-padding,1-padding);
        canv->GetPad(1)->SetLeftMargin(.11);
        canv->GetPad(1)->SetRightMargin(.05);
        //canv->GetPad(1)->SetBottomMargin(.3);
        canv->GetPad(1)->SetBottomMargin(.1);
        canv->GetPad(1)->Modified();

        canv->cd(1);
        gPad->Modified();
    }
}

// Basic settings for the histograms
void PlotMaker::setupHists(TH1D* data, TH1* ratio, TH1* ratio_1sig, TH1* ratio_2sig) {
    data->SetLineColor(1);
    data->SetLineWidth(2);
    data->SetMarkerStyle(20);
    data->SetMarkerSize(0.75);

    data->SetTitle(";;Events");
    data->GetYaxis()->SetTitleSize(0.05);
    data->GetYaxis()->SetTitleOffset(1.1);
    data->GetYaxis()->SetLabelSize(0.05);

    ratio->SetLineColor(1);     // same as h_data
    ratio->SetLineWidth(2);     // same as h_data
    ratio->SetMarkerStyle(20);  // same as h_data
    ratio->SetMarkerSize(0.75); // same as h_data
    //ratio->SetMarkerStyle(1);

    ratio_1sig->Reset();
    ratio_1sig->SetLineColor(1);     // same as h_data
    ratio_1sig->SetLineWidth(2);     // same as h_data
    ratio_1sig->SetMarkerStyle(20);  // same as h_data
    ratio_1sig->SetMarkerSize(0.75); // same as h_data
    ratio_1sig->SetMarkerStyle(1);
    ratio_1sig->SetFillColor(kGreen);
    ratio_1sig->SetLineColor(0);
    ratio_1sig->SetLineWidth(0);
    ratio_1sig->SetMarkerStyle(0);
    ratio_1sig->SetMarkerColor(0);
    ratio_1sig->SetMarkerSize(0.);

    ratio_2sig->Reset();
    ratio_2sig->SetLineColor(1);     // same as h_data
    ratio_2sig->SetLineWidth(2);     // same as h_data
    ratio_2sig->SetMarkerStyle(20);  // same as h_data
    ratio_2sig->SetMarkerSize(0.75); // same as h_data
    ratio_2sig->SetMarkerStyle(1);
    ratio_2sig->SetFillColor(kGreen+2);//kYellow
    ratio_2sig->SetLineColor(0);
    ratio_2sig->SetLineWidth(0);
    ratio_2sig->SetMarkerStyle(0);
    ratio_2sig->SetMarkerColor(0);
    ratio_2sig->SetMarkerSize(0.);

    // ratio_1sig is the histogram which gets drawn first, so only need to apply these settings once
    ratio_1sig->SetTitle(";channel and category");
    //ratio_1sig->GetYaxis()->SetNdivisions(50000+404);
    ratio_1sig->GetYaxis()->SetNdivisions(4,4,5);
    ratio_1sig->GetYaxis()->CenterTitle();
    ratio_1sig->GetYaxis()->SetLabelSize(0.1);
    ratio_1sig->GetYaxis()->SetTitle("Data/MC");
    ratio_1sig->GetYaxis()->SetTitleSize(0.1);
    ratio_1sig->GetYaxis()->SetTitleOffset(.45);
    ratio_1sig->GetXaxis()->SetLabelSize(0.1);
    ratio_1sig->GetXaxis()->SetTitleOffset(1.1);
    ratio_1sig->GetXaxis()->SetTitle(data->GetXaxis()->GetTitle());
    ratio_1sig->GetXaxis()->SetLabelSize(0.12);
    ratio_1sig->GetXaxis()->SetLabelOffset(0.02);
    ratio_1sig->GetXaxis()->SetTitleSize(0.13);
}

// Basic settings for the legends
void PlotMaker::setupLegends(TLegend *leg, TLegend *leg_ratio) {
    leg->SetFillColor(kWhite);
    leg->SetLineColor(kWhite);
    leg->SetShadowColor(kWhite);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);
    leg->SetNColumns(4);

    leg_ratio->SetFillColor(kWhite);
    leg_ratio->SetLineColor(kWhite);
    leg_ratio->SetShadowColor(kWhite);
    leg_ratio->SetTextFont(42);
    leg_ratio->SetTextSize(0.10);
}

#endif
/* PLOTMAKER */