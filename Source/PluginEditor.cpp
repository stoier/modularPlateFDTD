/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
PlateAudioProcessorEditor::PlateAudioProcessorEditor (PlateAudioProcessor& p)
    : AudioProcessorEditor (&p), audioProcessor (p),
 sig0SliderAttachment(audioProcessor.tree, "Frequency Independent Damping", sig0Slider), sig1SliderAttachment(audioProcessor.tree, "Frequency Dependent Damping", sig1Slider),
lengthXSliderAttachment(audioProcessor.tree, "Plate length X", lengthXSlider),
lengthYSliderAttachment(audioProcessor.tree, "Plate length Y", lengthYSlider),
excXSliderAttachment(audioProcessor.tree, "Excitation pos X", excXSlider),
excYSliderAttachment(audioProcessor.tree, "Excitation pos Y", excYSlider),
lisXSliderAttachment(audioProcessor.tree, "Listening pos X", lisXSlider),
lisYSliderAttachment(audioProcessor.tree, "Listening pos Y", lisYSlider),
thicknessSliderAttachment(audioProcessor.tree, "Plate thickness", thicknessSlider),
excFSliderAttachment(audioProcessor.tree, "Excitation force", excFSlider),
excTSliderAttachment(audioProcessor.tree, "Excitation time", excTSlider), vBSliderAttachment(audioProcessor.tree, "Bow velocity", vBSlider), fBSliderAttachment(audioProcessor.tree, "Bow force", fBSlider), fricSliderAttachment(audioProcessor.tree, "Friction", fricSlider), bAtt1SliderAttachment(audioProcessor.tree, "Bow attack 1", bAtt1Slider), bDec1SliderAttachment(audioProcessor.tree, "Bow decay 1", bDec1Slider), bSus1SliderAttachment(audioProcessor.tree, "Bow sustain 1", bSus1Slider), bRel1SliderAttachment(audioProcessor.tree, "Bow release 1", bRel1Slider), FBEnv1SliderAttachment(audioProcessor.tree, "Bow force env 1", FBEnv1Slider), vBEnv1SliderAttachment(audioProcessor.tree, "Bow velocity env 1", vBEnv1Slider), xPosLFOModSliderAttachment(audioProcessor.tree, "X Pos Mod Depth", xPosLFOModSlider), yPosLFOModSliderAttachment(audioProcessor.tree, "Y Pos Mod Depth", yPosLFOModSlider), lfoRateSliderAttachment(audioProcessor.tree, "LFO Rate", lfoRateSlider)
{

    
    setSliderAndLabelRotary(sig0Slider, sig0Label, "Freq Ind Damp");
    addAndMakeVisible(sig0Slider);
    addAndMakeVisible(sig0Label);
    
    setSliderAndLabelRotary(sig1Slider, sig1Label, "Freq Dep Damp");
    addAndMakeVisible(sig1Slider);
    addAndMakeVisible(sig1Label);
    
    
    setSliderAndLabelHorizontal(lengthXSlider, lengthXLabel, "Plate Length X");
    lengthXSlider.setTextValueSuffix(" m");
    addAndMakeVisible(lengthXSlider);
    addAndMakeVisible(lengthXLabel);
    
    setSliderAndLabelVertical(lengthYSlider, lengthYLabel, "Plate Length Y");
    lengthYSlider.setTextValueSuffix(" m");
    addAndMakeVisible(lengthYSlider);
    addAndMakeVisible(lengthYLabel);
    
    setSliderAndLabelHorizontal(excXSlider, excXLabel, "Excitation Pos X Ratio");
    addAndMakeVisible(excXSlider);
    addAndMakeVisible(excXLabel);
    
    setSliderAndLabelVertical(excYSlider, excYLabel, "Excitation Pos Y Ratio");
    addAndMakeVisible(excYSlider);
    addAndMakeVisible(excYLabel);
    
    setSliderAndLabelHorizontal(lisXSlider, lisXLabel, "Listening Pos X Ratio");
    addAndMakeVisible(lisXSlider);
    addAndMakeVisible(lisXLabel);
    
    setSliderAndLabelVertical(lisYSlider, lisYLabel, "Listening Pos Y Ratio");
    addAndMakeVisible(lisYSlider);
    addAndMakeVisible(lisYLabel);
    
    setSliderAndLabelRotary(thicknessSlider, thicknessLabel, "Thickness");
    thicknessSlider.setTextValueSuffix(" mm");
    addAndMakeVisible(thicknessSlider);
    addAndMakeVisible(thicknessLabel);
    
    setSliderAndLabelRotary(excFSlider, excFLabel, "Excitation Force");
    excFSlider.setTextValueSuffix(" N");
    addAndMakeVisible(excFSlider);
    addAndMakeVisible(excFLabel);
    
    setSliderAndLabelRotary(excTSlider, excTLabel, "Excitation Time");
    excTSlider.setTextValueSuffix(" ms");
    addAndMakeVisible(excTSlider);
    addAndMakeVisible(excTLabel);
    
    setSliderAndLabelRotary(vBSlider, vBLabel, "Bow Velocity");
    vBSlider.setTextValueSuffix(" m/s");
    addAndMakeVisible(vBSlider);
    addAndMakeVisible(vBLabel);
    
    setSliderAndLabelRotary(fBSlider, fBLabel, "Bow Force");
    excTSlider.setTextValueSuffix(" N");
    addAndMakeVisible(fBSlider);
    addAndMakeVisible(fBLabel);
    
    setSliderAndLabelRotary(fricSlider, fricLabel, "Friction");
    //fricSlider.setTextValueSuffix("");
    addAndMakeVisible(fricSlider);
    addAndMakeVisible(fricLabel);
    
    setSliderAndLabelVertical(bAtt1Slider, bAtt1Label, "Attack");
    addAndMakeVisible(bAtt1Slider);
    addAndMakeVisible(bAtt1Label);
    
    setSliderAndLabelVertical(bDec1Slider, bDec1Label, "Decay");
    addAndMakeVisible(bDec1Slider);
    addAndMakeVisible(bDec1Label);
    
    setSliderAndLabelVertical(bSus1Slider, bSus1Label, "Sustain");
    addAndMakeVisible(bSus1Slider);
    addAndMakeVisible(bSus1Label);
    
    setSliderAndLabelVertical(bRel1Slider, bRel1Label, "Release");
    addAndMakeVisible(bRel1Slider);
    addAndMakeVisible(bRel1Label);
    
    setSliderAndLabelRotary(FBEnv1Slider, FBEnv1Label, "Force Env Amount");
    addAndMakeVisible(FBEnv1Slider);
    addAndMakeVisible(FBEnv1Label);
    
    setSliderAndLabelRotary(vBEnv1Slider, vBEnv1Label, "Velocity Env Amount");
    addAndMakeVisible(vBEnv1Slider);
    addAndMakeVisible(vBEnv1Label);
    
    setSliderAndLabelRotary(lfoRateSlider, lfoRateLabel, "LFO Rate");
    addAndMakeVisible(lfoRateSlider);
    addAndMakeVisible(lfoRateLabel);
    
    setSliderAndLabelRotary(xPosLFOModSlider, xPosLFOModLabel, "X Pos Mod");
    addAndMakeVisible(xPosLFOModSlider);
    addAndMakeVisible(xPosLFOModLabel);
    
    setSliderAndLabelRotary(yPosLFOModSlider, yPosLFOModLabel, "Y Pos Mod");
    addAndMakeVisible(yPosLFOModSlider);
    addAndMakeVisible(yPosLFOModLabel);
    
    
    hitButton.setColour(juce::TextButton::ColourIds::buttonOnColourId, juce::Colours::purple);
    addAndMakeVisible(hitButton);
    hitButton.onClick = [this] {hitButtonClicked();};
    
    bowExcButton.setColour(juce::TextButton::ColourIds::buttonOnColourId, juce::Colours::purple);
    bowExcButton.setRadioGroupId(1);
    bowExcButton.setToggleable(true);
    bowExcButton.setClickingTogglesState(true);
    //bowExcButton.setToggleState(true);
    addAndMakeVisible(bowExcButton);
    bowExcButton.onStateChange = [this]
    {
        if (bowExcButton.getToggleState() == true)
        {
            resized();
            audioProcessor.excTypeId = 1;
        }
            
    };
    
    malletExcButton.setColour(juce::TextButton::ColourIds::buttonOnColourId, juce::Colours::purple);
    malletExcButton.setRadioGroupId(1);
    malletExcButton.setToggleable(true);
    malletExcButton.setClickingTogglesState(true);
    addAndMakeVisible(malletExcButton);
    malletExcButton.onStateChange = [this]
    {
        if (malletExcButton.getToggleState() == true)
        {
            resized();
            audioProcessor.excTypeId = 2;
        }
            
    };
    
    /*
    linkFBvB.setColour(juce::TextButton::ColourIds::buttonOnColourId, juce::Colours::purple);
    linkFBvB.setToggleable(true);
    linkFBvB.setClickingTogglesState(true);
    addAndMakeVisible(linkFBvB);
    linkFBvB.onStateChange = [this]
    {
        if (linkFBvB.getToggleState() == true)
        {
           
        }
            
    };*/

    startBowButton.setColour(juce::TextButton::ColourIds::buttonOnColourId, juce::Colours::purple);
    startBowButton.setToggleable(true);
    startBowButton.setClickingTogglesState(true);
    addAndMakeVisible(startBowButton);
    startBowButton.onStateChange = [this]
    {
        if (startBowButton.getToggleState() == true)
        {
            audioProcessor.bowEnd = false;
            audioProcessor.bowStart = true;
        }
        else
        {
            audioProcessor.bowEnd = true;
            audioProcessor.bowStart = false;
        }
            
    };
    
    plateMaterialMenu.addItem("Brass", 1);
    plateMaterialMenu.addItem("Bronze", 2);
    plateMaterialMenu.addItem("Iron", 3);
    plateMaterialMenu.addItem("Aluminium", 4);
    plateMaterialMenu.addItem("Gold", 5);
    plateMaterialMenu.addItem("Silver", 6);
    plateMaterialMenu.addItem("Copper", 7);
    plateMaterialMenu.setSelectedId (audioProcessor.plateMaterialId);
    plateMaterialMenu.setSelectedId(1);
    plateMaterialMenu.onChange = [this]
    {
        switch (plateMaterialMenu.getSelectedId())
        {
            case 1: audioProcessor.plateMaterialId = 1; break;
            case 2: audioProcessor.plateMaterialId = 2; break;
            case 3: audioProcessor.plateMaterialId = 3; break;
            case 4: audioProcessor.plateMaterialId = 4; break;
            case 5: audioProcessor.plateMaterialId = 5; break;
            case 6: audioProcessor.plateMaterialId = 6; break;
            case 7: audioProcessor.plateMaterialId = 7; break;
            default: break;
        };
    };
    addAndMakeVisible(plateMaterialMenu);
    plateMaterialMenuLabel.setText("Plate Material:", juce::dontSendNotification);
    plateMaterialMenuLabel.setFont(18.0f);
    plateMaterialMenuLabel.setJustificationType(juce::Justification::left);
    plateMaterialMenuLabel.attachToComponent(&plateMaterialMenu , false);
    
    setResizable(true, true);
    setResizeLimits(400, 250, 1600, 1000);
    getConstrainer() -> setFixedAspectRatio(1.2);
    
    
    setSize (800, 500);
}

PlateAudioProcessorEditor::~PlateAudioProcessorEditor()
{
}

void PlateAudioProcessorEditor::hitButtonClicked()
{
    audioProcessor.hit = true;
}


void PlateAudioProcessorEditor::setSliderAndLabelRotary(juce::Slider& slider, juce::Label& label, const juce::String &Text)
{
    slider.setSliderStyle(juce::Slider::SliderStyle::RotaryVerticalDrag);
    slider.setColour(0x1001311, juce::Colours::purple);
    slider.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 60, 20);
    label.setText(Text, juce::dontSendNotification);
    label.setFont(15.0f);
    label.setJustificationType(juce::Justification::centred);
    label.attachToComponent(&slider , false);
}


void PlateAudioProcessorEditor::setSliderAndLabelHorizontal(juce::Slider& slider, juce::Label& label, const juce::String &Text)
{
    slider.setSliderStyle(juce::Slider::SliderStyle::LinearHorizontal);
    slider.setColour(0x1001310, juce::Colours::purple);
    slider.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 60, 20);
    label.setText(Text, juce::dontSendNotification);
    label.setFont(15.0f);
    label.setJustificationType(juce::Justification::centred);
    label.attachToComponent(&slider , false);
}


void PlateAudioProcessorEditor::setSliderAndLabelVertical(juce::Slider& slider, juce::Label& label, const juce::String &Text)
{
    slider.setSliderStyle(juce::Slider::SliderStyle::LinearVertical);
    slider.setColour(0x1001310, juce::Colours::purple);
    slider.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 60, 20);
    label.setText(Text, juce::dontSendNotification);
    label.setFont(15.0f);
    label.setJustificationType(juce::Justification::centred);
    label.attachToComponent(&slider , false);
}


//==============================================================================
void PlateAudioProcessorEditor::paint (juce::Graphics& g)
{
    g.fillAll (getLookAndFeel().findColour (juce::ResizableWindow::backgroundColourId));

    g.fillAll (juce::Colours::black);
    g.setColour (juce::Colours::white);
    g.setFont(20.f);
}

void PlateAudioProcessorEditor::resized()
{
    /*
    auto excitationArea = getLocalBounds();
    auto plateArea = excitationArea.removeFromTop(excitationArea.getHeight()*0.4);
    plateArea.removeFromTop(plateArea.getHeight()*0.15);
    auto dampingArea0 = plateArea.removeFromLeft(plateArea.getWidth()*0.35);
    dampingArea0.removeFromBottom(dampingArea0.getHeight()*0.4);
    auto dampingArea1 = dampingArea0.removeFromLeft(dampingArea0.getWidth()*0.5);
    sig0Slider.setBounds(dampingArea0);
    sig1Slider.setBounds(dampingArea1);
    auto ThicknessArea =plateArea.removeFromLeft(plateArea.getWidth()*0.25);
    thicknessSlider.setBounds(ThicknessArea.removeFromTop(ThicknessArea.getHeight()*0.6));
    plateArea.removeFromBottom(plateArea.getHeight()*0.3);
    lengthXSlider.setBounds(plateArea.removeFromLeft(plateArea.getWidth()*0.5));
    lengthYSlider.setBounds(plateArea);
    auto excitationArea2 = excitationArea.removeFromTop(excitationArea.getHeight()*0.5);
    excitationArea.removeFromBottom(excitationArea.getHeight()*0.2);
    excXSlider.setBounds(excitationArea.removeFromLeft(excitationArea.getWidth()*0.25));
    excYSlider.setBounds(excitationArea.removeFromLeft(excitationArea.getWidth()*0.33));
    lisXSlider.setBounds(excitationArea.removeFromLeft(excitationArea.getWidth()*0.5));
    lisYSlider.setBounds(excitationArea);
    excitationArea2.removeFromBottom(excitationArea2.getHeight()*0.3);
    excFSlider.setBounds(excitationArea2.removeFromLeft(excitationArea2.getWidth()*0.25));
    excTSlider.setBounds(excitationArea2.removeFromLeft(excitationArea2.getWidth()*0.33));
    hitButton.setBounds(excitationArea2.removeFromLeft(excitationArea2.getWidth()*0.33));
    excitationArea2.removeFromTop(excitationArea2.getHeight()*0.2);
    excitationArea2.removeFromLeft(excitationArea2.getWidth()*0.2);
    excitationArea2.removeFromRight(excitationArea2.getWidth()*0.1);
    plateMaterialMenu.setBounds(excitationArea2);
    */
    auto excitationArea = getLocalBounds();
    auto plateArea = excitationArea.removeFromRight(excitationArea.getWidth()*0.3);
    auto dampingArea = plateArea.removeFromTop(plateArea.getHeight()*0.25);
    auto dampingAreaHeigth = dampingArea.getHeight();
    dampingArea.removeFromTop(dampingAreaHeigth*0.3);
    auto rotaryWidth = dampingArea.getWidth()*0.33;
    sig0Slider.setBounds(dampingArea.removeFromLeft(rotaryWidth));
    sig1Slider.setBounds(dampingArea.removeFromLeft(rotaryWidth));
    auto thicknessArea = dampingArea.removeFromRight(dampingArea.getWidth());
    auto plateXYArea = plateArea.removeFromTop(dampingAreaHeigth*1.2);
    plateXYArea.removeFromTop(dampingAreaHeigth*0.3);
    thicknessSlider.setBounds(thicknessArea.removeFromTop(dampingAreaHeigth));
    lengthXSlider.setBounds(plateXYArea.removeFromLeft(plateXYArea.getWidth()*0.5));
    lengthYSlider.setBounds(plateXYArea);
    auto listeningArea = plateArea.removeFromTop(dampingAreaHeigth*1.2);
    listeningArea.removeFromTop(dampingAreaHeigth*0.3);
    lisXSlider.setBounds(listeningArea.removeFromLeft(listeningArea.getWidth()*0.5));
    lisYSlider.setBounds(listeningArea);
    plateArea.removeFromTop(dampingAreaHeigth*0.3);
    plateMaterialMenu.setBounds(plateArea);
    
    auto excitationTypeArea = excitationArea.removeFromTop(excitationArea.getHeight()*0.1);
    bowExcButton.setBounds(excitationTypeArea.removeFromLeft(excitationArea.getWidth()*0.5));
    malletExcButton.setBounds(excitationTypeArea);
    if (bowExcButton.getToggleState() == true)
    {
        hitButton.setVisible(false);
        excTSlider.setVisible(false);
        excFSlider.setVisible(false);
        vBSlider.setVisible(true);
        fBSlider.setVisible(true);
        startBowButton.setVisible(true);
        fricSlider.setVisible(true);
        bAtt1Slider.setVisible(true);
        bDec1Slider.setVisible(true);
        bSus1Slider.setVisible(true);
        bRel1Slider.setVisible(true);
        FBEnv1Slider.setVisible(true);
        vBEnv1Slider.setVisible(true);
        lfoRateSlider.setVisible(true);
        xPosLFOModSlider.setVisible(true);
        yPosLFOModSlider.setVisible(true);
        
        auto ffvBArea = excitationArea.removeFromTop(dampingAreaHeigth);
        ffvBArea.removeFromTop(dampingAreaHeigth*0.3);
        vBSlider.setBounds(ffvBArea.removeFromLeft(excitationArea.getWidth()*0.33));
        fBSlider.setBounds(ffvBArea.removeFromLeft(excitationArea.getWidth()*0.5));
        fricSlider.setBounds(ffvBArea);
        
        auto ADSR11Area = excitationArea.removeFromBottom(excitationArea.getHeight()*0.5);
        ADSR11Area.removeFromTop(dampingAreaHeigth*0.3);
        auto LFOArea = ADSR11Area.removeFromRight(ADSR11Area.getWidth()*0.5);
        auto envSliderWidth = 0.25 * ADSR11Area.getWidth();
        bAtt1Slider.setBounds(ADSR11Area.removeFromLeft(envSliderWidth));
        bDec1Slider.setBounds(ADSR11Area.removeFromLeft(envSliderWidth));
        bSus1Slider.setBounds(ADSR11Area.removeFromLeft(envSliderWidth));
        bRel1Slider.setBounds(ADSR11Area.removeFromLeft(envSliderWidth));
        auto lfoSliderWidth = 0.33 * LFOArea.getWidth();
        lfoRateSlider.setBounds(LFOArea.removeFromLeft(lfoSliderWidth));
        xPosLFOModSlider.setBounds(LFOArea.removeFromLeft(lfoSliderWidth));
        yPosLFOModSlider.setBounds(LFOArea.removeFromLeft(lfoSliderWidth));
        //FBEnv1Slider.setBounds(ADSR12Area.removeFromLeft(ADSR12Area.getWidth()*0.5));
        //vBEnv1Slider.setBounds(ADSR12Area);
        excitationArea.removeFromTop(dampingAreaHeigth*0.3);
        auto startBowArea = excitationArea.removeFromRight(excitationArea.getWidth()*0.33);
        startBowButton.setBounds(startBowArea);
        excXSlider.setBounds(excitationArea.removeFromLeft(excitationArea.getWidth()*0.5));
        excYSlider.setBounds(excitationArea);
    }
    
    if (malletExcButton.getToggleState() == true)
    {
        vBSlider.setVisible(false);
        fBSlider.setVisible(false);
        fricSlider.setVisible(false);
        bAtt1Slider.setVisible(false);
        bDec1Slider.setVisible(false);
        bSus1Slider.setVisible(false);
        bRel1Slider.setVisible(false);
        FBEnv1Slider.setVisible(false);
        vBEnv1Slider.setVisible(false);
        lfoRateSlider.setVisible(false);
        xPosLFOModSlider.setVisible(false);
        yPosLFOModSlider.setVisible(false);
        startBowButton.setVisible(false);
        hitButton.setVisible(true);
        excTSlider.setVisible(true);
        excFSlider.setVisible(true);
        
        auto excFTArea = excitationArea.removeFromTop(dampingAreaHeigth);
        excFTArea.removeFromTop(dampingAreaHeigth*0.3);
        excFSlider.setBounds(excFTArea.removeFromLeft(excitationArea.getWidth()*0.5));
        excTSlider.setBounds(excFTArea);
        auto hitArea = excitationArea.removeFromTop(dampingAreaHeigth);
        hitArea.removeFromTop(dampingAreaHeigth*0.3);
        hitButton.setBounds(hitArea);
        excitationArea.removeFromTop(dampingAreaHeigth*0.3);
        excXSlider.setBounds(excitationArea.removeFromLeft(excitationArea.getWidth()*0.5));
        excYSlider.setBounds(excitationArea);
    }

}
