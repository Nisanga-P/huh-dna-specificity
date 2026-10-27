import torch.nn as nn


class BinaryClassiferNN(nn.Module):
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

        layers.append(nn.Linear(current_dim, 1))
        self.network = nn.Sequential(*layers)

    def forward(self, x):
        return self.network(x)