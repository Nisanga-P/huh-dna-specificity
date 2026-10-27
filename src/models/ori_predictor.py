import torch.nn as nn


class MultiHeadFourClassClassiferNN(nn.Module):
    def __init__(
        self,
        input_dim,
        hidden_dims,
        activation=nn.ReLU,
        dropout_rate=0.3,
    ):
        super().__init__()

        if isinstance(hidden_dims, int):
            hidden_dims = [hidden_dims]

        layers = []
        current_dim = input_dim

        for hidden_dim in hidden_dims:
            layers.append(nn.Linear(current_dim, hidden_dim))
            layers.append(activation())
            if dropout_rate > 0:
                layers.append(nn.Dropout(dropout_rate))
            current_dim = hidden_dim

        self.shared_backbone = nn.Sequential(*layers)
        self.heads = nn.ModuleList(
            [nn.Linear(current_dim, 4) for _ in range(9)]
        )

    def forward(self, x):
        x = self.shared_backbone(x)
        return [head(x) for head in self.heads]