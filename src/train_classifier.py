import torch
import torch.nn as nn

from data_setup import train_dataloader, X_val_rep, X_val_relax, y_val_rep, y_val_relax
from models import BinaryClassiferNN, BinaryClassifierTwoHeadNN


SEED = 3242
EPOCHS = 100
LEARNING_RATE = 0.001
MODEL_TYPE = "single_input"
DNA_DIM = 9

SINGLE_INPUT_DIM = 34
SINGLE_HIDDEN_DIMS = [48, 32]

TWO_HEAD_SDBM_DIM = 118
TWO_HEAD_DNA_DIM = 9
TWO_HEAD_HIDDEN_DIM = 32
TWO_HEAD_DIM = 24
DROPOUT_P = 0.2


torch.manual_seed(SEED)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

if MODEL_TYPE == "single_input":
    protein_dna_bind_classifier = BinaryClassiferNN(
        input_dim=SINGLE_INPUT_DIM,
        hidden_dims=SINGLE_HIDDEN_DIMS,
    ).to(device)
elif MODEL_TYPE == "two_head":
    protein_dna_bind_classifier = BinaryClassifierTwoHeadNN(
        sdbm_dim=TWO_HEAD_SDBM_DIM,
        dna_dim=TWO_HEAD_DNA_DIM,
        hidden_dim1=TWO_HEAD_HIDDEN_DIM,
        head_dim=TWO_HEAD_DIM,
        dropout_p=DROPOUT_P,
    ).to(device)
else:
    raise ValueError(
        f"Unknown MODEL_TYPE: {MODEL_TYPE!r}. "
        "Choose 'single_input' or 'two_head'."
    )

loss_fn = nn.BCEWithLogitsLoss()

optimiser = torch.optim.Adam(
    protein_dna_bind_classifier.parameters(),
    lr=LEARNING_RATE,
)


def accuracy_fn(y_true, y_pred):
    correct_pred = torch.eq(y_true, y_pred).sum().item()
    return (correct_pred / len(y_pred)) * 100


x_test_dataset = torch.cat((X_val_rep, X_val_relax), dim=0)
y_test_dataset = torch.cat((y_val_rep, y_val_relax), dim=0)

acc = []
acc_val = []
loss_train = []
loss_val = []

for epoch in range(EPOCHS):
    protein_dna_bind_classifier.train()

    running_loss = 0.0
    running_acc = 0.0

    for inputs, targets in train_dataloader:
        inputs = inputs.to(device, non_blocking=True)
        targets = targets.to(device, non_blocking=True).float().unsqueeze(dim=1)

        if MODEL_TYPE == "single_input":
            logits = protein_dna_bind_classifier(inputs)
        else:
            sdbm_input = inputs[:, :-DNA_DIM]
            dna_input = inputs[:, -DNA_DIM:]
            logits = protein_dna_bind_classifier(sdbm_input, dna_input)

        predictions = torch.round(torch.sigmoid(logits))
        accuracy = accuracy_fn(targets, predictions)
        loss = loss_fn(logits, targets)

        optimiser.zero_grad()
        loss.backward()
        optimiser.step()

        running_loss += loss.item()
        running_acc += accuracy

    epoch_loss = running_loss / len(train_dataloader)
    epoch_acc = running_acc / len(train_dataloader)

    protein_dna_bind_classifier.eval()

    with torch.inference_mode():
        test_inputs = x_test_dataset.to(device, non_blocking=True)
        test_targets = (
            y_test_dataset.to(device, non_blocking=True)
            .float()
            .unsqueeze(dim=1)
        )

        if MODEL_TYPE == "single_input":
            test_logits = protein_dna_bind_classifier(test_inputs)
        else:
            test_sdbm_input = test_inputs[:, :-DNA_DIM]
            test_dna_input = test_inputs[:, -DNA_DIM:]
            test_logits = protein_dna_bind_classifier(
                test_sdbm_input,
                test_dna_input,
            )

        test_predictions = torch.round(torch.sigmoid(test_logits))
        test_loss = loss_fn(test_logits, test_targets)
        test_accuracy = accuracy_fn(test_targets, test_predictions)

    acc.append(epoch_acc)
    acc_val.append(test_accuracy)
    loss_train.append(epoch_loss)
    loss_val.append(test_loss.item())

    print(
        f"Epoch: {epoch + 1:02d}/{EPOCHS} | "
        f"Training Loss: {epoch_loss:.4f} | "
        f"Training Accuracy: {epoch_acc:.4f} | "
        f"Validation Loss: {test_loss.item():.4f} | "
        f"Validation Accuracy: {test_accuracy:.4f}"
    )