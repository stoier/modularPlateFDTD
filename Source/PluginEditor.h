/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "PluginProcessor.h"
#include "ThinPlate.h"


//==============================================================================
/**
*/
class PlateAudioProcessorEditor  : public juce::AudioProcessorEditor
{
public:
    PlateAudioProcessorEditor (PlateAudioProcessor&);
    ~PlateAudioProcessorEditor() override;

    //==============================================================================
    void paint (juce::Graphics&) override;
    void resized() override;
    
    void setSliderAndLabelRotary(juce::Slider& slider, juce::Label& label, const juce::String &Text);

    void setSliderAndLabelHorizontal(juce::Slider& slider, juce::Label& label, const juce::String &Text);

    void setSliderAndLabelVertical(juce::Slider& slider, juce::Label& label, const juce::String &Text);

    void hitButtonClicked();
    
    void bowExcButtonClicked();
    
    void malletExcButtonClicked();
    
private:

    PlateAudioProcessor& audioProcessor;
    
    using APVTS = juce::AudioProcessorValueTreeState;
    using sliderAttachment = APVTS::SliderAttachment;
    
    juce::Slider sig0Slider, sig1Slider, lengthXSlider, lengthYSlider, excXSlider, excYSlider, lisXSlider, lisYSlider, excFSlider, excTSlider, thicknessSlider, vBSlider, fBSlider, fricSlider, bAtt1Slider, bDec1Slider, bSus1Slider, bRel1Slider, FBEnv1Slider, vBEnv1Slider, xPosLFOModSlider, yPosLFOModSlider, lfoRateSlider;
    
    juce::Label sig0Label, sig1Label, lengthXLabel, lengthYLabel, excXLabel, excYLabel, lisXLabel, lisYLabel, excFLabel, excTLabel, thicknessLabel, plateMaterialMenuLabel, vBLabel, fBLabel, fricLabel, bAtt1Label, bDec1Label, bSus1Label, bRel1Label, FBEnv1Label, vBEnv1Label, xPosLFOModLabel, yPosLFOModLabel, lfoRateLabel;
    
    juce::ComboBox plateMaterialMenu;
    
    sliderAttachment sig0SliderAttachment, sig1SliderAttachment, lengthXSliderAttachment, lengthYSliderAttachment, excXSliderAttachment, excYSliderAttachment, lisXSliderAttachment, lisYSliderAttachment, thicknessSliderAttachment, excFSliderAttachment, excTSliderAttachment, vBSliderAttachment, fBSliderAttachment, fricSliderAttachment, bAtt1SliderAttachment, bDec1SliderAttachment, bSus1SliderAttachment, bRel1SliderAttachment, FBEnv1SliderAttachment, vBEnv1SliderAttachment, xPosLFOModSliderAttachment, yPosLFOModSliderAttachment, lfoRateSliderAttachment;
    
    
    juce::TextButton hitButton{"Hit plate"}, bowExcButton{"Bow"}, malletExcButton{"Mallet"}, linkFBvB{"Link"}, startBowButton{"Start bowing"};


    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (PlateAudioProcessorEditor)
};

