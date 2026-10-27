import torch
import torch.nn as nn


class BinaryClassifierTwoHeadNN(nn.Module):
    def __init__(
        self,
        sdbm_dim,
        dna_dim=9,
        hidden_dim1=32,
        head_dim=24,
        dropout_p=0.2,
    ):
        super().__init__()

        self.sdbm_fc = nn.Linear(sdbm_dim, head_dim)
        self.dna_fc = nn.Linear(dna_dim, head_dim)

        self.relu = nn.ReLU()
        self.dropout = nn.Dropout(p=dropout_p)

        self.fc1 = nn.Linear(head_dim * 2, hidden_dim1)
        self.fc2 = nn.Linear(hidden_dim1, 1)

    def forward(self, sdbm_input, dna_input):
        sdbm = self.dropout(self.relu(self.sdbm_fc(sdbm_input)))
        dna = self.dropout(self.relu(self.dna_fc(dna_input)))

        x = torch.cat((sdbm, dna), dim=1)
        x = self.dropout(self.relu(self.fc1(x)))

        return self.fc2(x)