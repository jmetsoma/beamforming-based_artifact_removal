# beamforming-based_artifact_removal
These Matlab codes allow you to design spatial filter for artifact removal from EEG or MEG

List of functions

BF_type_cleaningFilt.m - the basic function for using beamforming- based EEG cleaning, when you know the topographies of the artifact subspace

SOUND_veryfast.m - SOUND cleaning algorithm implemented via beamforming, which requires no looping over channels, making it fast

simple_wiener_veryfast.m - used by SOUND function ( or can be used separately), data driven version of SOUND

FilterOutGivenPatternsSound_BF.m - beamformin based cleaning with given topography, when the input data were cleaned by SOUND and topography obtained before SOUND

SOUND_fast_SSP.m - joint SOUND and SSP-SIR algorithm, for simultaneous usage, where SSP is embedded in SOUND

removeShamBF.m - removing from an interesting EEG condition (active TMS or some other condition) the signal produced by sources which generated sham EEG (or some other condition) 


