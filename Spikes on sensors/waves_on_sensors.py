import mne
from mne.preprocessing import ICA

Data = mne.io.read_raw_ctf('/home/sasha/Data/Wave_prior/Data/MEG/Kaco/Kaco/kaco_spont01.ds')

Data.load_data()

Data.info

picks = mne.pick_types(Data.info, meg='mag', eeg=False, eog=False,
                       stim=False, ref_meg = False, exclude='bads')
Data.filter(1, 40., fir_design='firwin', picks = picks)

n_components = 68
method = 'fastica'  
decim = 3  
random_state = 23
ica = ICA(n_components=n_components, method=method, random_state=random_state)
ica.fit(Data, picks=picks, decim=decim)

ica.plot_components(picks=range(40), inst=Data)

ica.plot_sources(Data)

ica.apply(Data, exclude = [0,3,12,19])

Data.save('/home/sasha/Data/Wave_prior/Data/MEG/Kaco/Kaco/kaco_spont01_py.fif')