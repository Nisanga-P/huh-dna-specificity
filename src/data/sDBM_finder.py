import torch
import torch.nn.functional as F
from esm.models.esmc import EsmcForMaskedLM, EsmcTokenizer

class sDBM_finder:

    '''
    Constructor for 'sDBM_finder' class.

    Choose model (default 300M variant of ESM-C), and specify reference sDBM to match against (PCV/WDV).
    
    '''
    def __init__(self, model_name = "ESMC-300M"):

        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

        model_str = "biohub/" + model_name

        self.model = EsmcForMaskedLM.from_pretrained(
                    model_str,
                    device=self.device,
                ).eval()

        self.tokeniser = EsmcTokenizer()


    # HELPER METHOD!
    def _get_protein_embedding(self, amino_acid_seq):

        inputs = self.tokeniser([amino_acid_seq], return_tensors="pt", padding=True)
        inputs = {name: value.to(self.model.device) for name, value in inputs.items()}

        with torch.inference_mode():
            result = self.model(**inputs, output_hidden_states=True)

        token_embeddings = result.hidden_states[-1]    
        residue_embeddings = token_embeddings[0, 1:-1] # exlude start and end token

        return residue_embeddings

    '''

    Obtains the sDBM embedding given known sDBM position in sequence (by using specified ESM-C model).
    
    Args:
     - amino_acid_seq: amino acid sequence of Rep protein to be used as reference as a str (with sDBM position in sequence already known)
     - sDBM_start_pos: starting position of sDBM in the provided amino acid sequence (indexing beginning at 1)

    Returns:
        reference sDBM embedding tensor 
    
    '''

    def get_reference_sDBM_embedding(self, amino_acid_seq, sDBM_start_pos):

        # fix indexing convention from 'start at 1' -> 'start at 0'
        sDBM_start_pos -= 1
        sDBM_end_pos = sDBM_start_pos + 9

        ref_protein_embed = self._get_protein_embedding(amino_acid_seq)

        return ref_protein_embed[sDBM_start_pos : sDBM_end_pos, :]


    '''
    
    Obtains sDBM for a Rep protein.

    Args:
     - amino_acid_sequence: sequence of Rep protein 
     - ref_sDBM_embedding: reference sDBM embedding used for match (output of 'get_reference_sDBM_embedding', likely PCV or WDV)
     - search_frame: range of positions in amino acid sequence that is searched over (default from 30 to 150) -> indexed from 1
    
    Return:
     tuple of (predicted sDBM sequence, start and end positions, confidence)

    
    '''

    def find_sDBM(self, amino_acid_sequence, ref_sDBM_embedding, search_frame = (30, 150)):

        protein_embed = self._get_protein_embedding(amino_acid_sequence)

        length, _ = protein_embed.shape

        start_range, end_range = search_frame

        # fix indexing convention from 'start at 1' -> 'start at 0'
        start_range -= 1
        end_range -= 1


        # amend start_range and end_range to actually acceptable values in the sequence
        if start_range < 0 or end_range > (length - 8):
            start_range = 0
            end_range = length - 8

        # mean of cosine similarity cannot go below -1
        max_mean = -1
        next_max_mean = -1

        start_seq = 0
        end_seq = 0 

        # iterate over each 9-length block within the specified window of embeddig, computing the mean cosine similarity at each point
        # max mean cos similarity used as the position for sDBM prediction
        # second max also recorded for purpose of a confidence metric
        for i in range(start_range, end_range):
            match_block = protein_embed[i:i+9, :]

            sim = F.cosine_similarity(match_block, ref_sDBM_embedding, dim = 1)

            mean = torch.mean(sim)

            if mean > max_mean:
                # reverse indexing convention from 'start at 0' -> 'start at 1'
                start_seq = i + 1
                end_seq = i + 9
                next_max_mean = max_mean
                max_mean = mean
            elif mean > next_max_mean:
                next_max_mean = mean

        # measure of confidence of the prediction
        diff_perc = (max_mean - next_max_mean) * 100

        sDBM_pred = str(amino_acid_sequence[start_seq - 1 : end_seq])

        sDBM_pos = (start_seq, end_seq)

        return (sDBM_pred, sDBM_pos, round(float(diff_perc), 1))





        





    

    
