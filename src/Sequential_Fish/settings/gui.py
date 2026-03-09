"""
Submodule handling data opening and merging for viewer
"""
from types import NoneType, UnionType
from PyQt5.QtWidgets import (
    QPushButton, QHBoxLayout, QLabel,
    QDialog, QFormLayout, QLineEdit, QSpinBox, QDoubleSpinBox,QCheckBox, QWidget, QPlainTextEdit
)
from typing import Dict, Any,get_args, get_origin, Type
from pandas import annotations
from pydantic import BaseModel


class ParametersModifier(QDialog):
    def __init__(self, data_model : Type[BaseModel], parent: QWidget | None = None, **parameters: Any):
        super().__init__(parent)
        self.setWindowTitle("Parameter Modification")
        print(parameters)
        self._data_model = data_model
        self._param_types = data_model.model_fields
        self._default_value = parameters
        self._widgets: Dict[str, Any] = {}

        form = QFormLayout()
       # Dynamically create widgets with default.types
        for name, default in self._param_types.items():
            att_type = default.annotation
            type_origin = get_origin(att_type)
            att_args = []

            widget = None
            add_none_option = False
            if type_origin is UnionType :
                allowed_types = list(get_args(default.annotation))
                if NoneType in allowed_types :
                    add_none_option = True
                    allowed_types.remove(NoneType)
                assert len(allowed_types) == 1, allowed_types
                att_type = allowed_types[0]
                type_origin = get_origin(att_type)
            
            if type_origin is None : #Simple types
                pass
            else : #types with args (tuples, list...)
                att_args = get_args(att_type)
                att_type = type_origin

            if att_type is int :
                sb = QSpinBox()
                sb.setRange(-1_000_000, 1_000_000)
                value = self._default_value[name]
                sb.setValue(value if not value is None else 0)
                widget = sb

            # Float self._default_value[name]
            elif att_type is float:
                sb = QDoubleSpinBox()
                sb.setRange(-1e6, 1e6)
                sb.setDecimals(6)
                value = self._default_value[name]
                sb.setValue(value if not value is None else 0)
                widget = sb

            # String self._default_value[name]
            elif att_type is str:
                le = QLineEdit()
                value = self._default_value[name]
                le.setText(value if not value is None else "")
                widget = le

            # Bool self._default_value[name]
            elif att_type is bool :
                cb = QCheckBox()
                cb.setChecked(self._default_value[name])
                widget = cb

            # Tuple of ints
            elif att_type is tuple :
                if all(x is int for x in att_args):
                    spinboxes = []
                    sub_layout = QHBoxLayout()
                    for val in self._default_value[name]:
                        sb = QSpinBox()
                        sb.setRange(-1_000_000, 1_000_000)
                        sb.setValue(val if not val is None else 0)
                        sub_layout.addWidget(sb)
                        spinboxes.append(sb)
                    widget = sub_layout
                    self._widgets[name + "__tuple__"] = spinboxes
            
                elif all(x is str for x in att_args):
                    line_edits = []
                    sub_layout = QHBoxLayout()
                    for val in self._default_value[name]:
                        le = QLineEdit()
                        le.setText(val)
                        sub_layout.addWidget(le)
                        line_edits.append(le)
                    widget = sub_layout
                    self._widgets[name + "__tuple__"] = line_edits
                else :
                    raise NotImplementedError()

            # List of simple types
            elif att_type is list :
                le = QLineEdit()
                # join list to comma-separated string
                le.setText(
                    ','.join(str(x) if not x is None else '' for x in self._default_value[name]) if not self._default_value[name] is None else ""
                )
                widget = le
            
            # Dict[str, simple]
            elif att_type is dict :
                pt = QPlainTextEdit()
                # key:value per line
                lines = [f"{k}:{v}" for k, v in self._default_value[name].items()]
                pt.setPlainText("\n".join(lines)); pt.setFixedHeight(80)
                widget = pt

            # Fallback: show as string
            else:
                raise NotImplementedError(f"not gui implementation for type {name}; {att_type} --> {type_origin} from {default.annotation}. {att_args}")

            form.addRow(QLabel(name.replace('_', ' ').capitalize()), widget)
            if add_none_option :
                sub_layout = QHBoxLayout()
                sub_layout.addWidget(widget)
                cb = QCheckBox()
                cb.setChecked(self._default_value[name] is None)
                sub_layout.addWidget(cb)
                self._widgets[name+"__tuple__"] = [widget, cb]

            self._widgets[name] = widget

        # Buttons
        btn_layout = QHBoxLayout()
        ok_btn = QPushButton("OK")
        ok_btn.clicked.connect(self.accept)
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        btn_layout.addWidget(ok_btn)
        btn_layout.addWidget(cancel_btn)

        form.addRow(btn_layout)
        self.setLayout(form)

    def get_parameters(self) -> Dict[str, Any]:
        """Retrieve parameters as a dict, converting widget values to Python types."""
        result: Dict[str, Any] = {}
        print(self._widgets)
        for name, default in self._param_types.items():

            if not name in self._widgets.keys() :
                name += "__tuple__"

            widget = self._widgets[name]
            att_type = default.annotation
            type_origin = get_origin(att_type)
            att_args = []
            add_none_option = False

            if type_origin is None : #Simple types
                pass
            else : #types with args (tuples, list...)
                att_args = get_args(att_type)
                att_type = type_origin

            if type_origin is UnionType :
                allowed_types = list(get_args(default.annotation))
                if NoneType in allowed_types :
                    add_none_option = True
                    allowed_types.remove(NoneType)
                assert len(allowed_types) == 1, allowed_types
                att_type = allowed_types[0]
                type_origin = get_origin(att_type)


            # int
            if add_none_option and widget[1].isChecked() :
                pass
            if att_type is int :
                result[name] = widget.value()
            # float
            elif att_type is float:
                result[name] = widget.value()
            # str
            elif att_type is str:
                result[name] = widget.text()
            # bool
            elif att_type is bool:
                result[name] = widget.isChecked()
            # tuple of ints
            elif att_type is tuple :
                if all(x is int for x in att_args):
                    spinboxes = self._widgets[name + "__tuple__"]
                    result[name] = tuple(sb.value() for sb in spinboxes)
                elif all(x is str for x in att_args):
                    line_edits = self._widgets[name + "__tuple__"]
                    result[name] = tuple(le.text() for le in line_edits)
                else : 
                    raise NotImplementedError()
            
            # List
            elif att_type is list:
                text = widget.text(); items = [i.strip() for i in text.split(',')]
                conv = type(self._default_value[name][0])
                result_items = []
                if text != "" :
                    for i in items :
                        try :
                            result_items.append(conv(i) if i.strip() != "" else None)
                        except ValueError :
                            print(f"Incorrect value set for {name} : expected {conv} got {i}.\nParameter restaured to {self._default_value[name]}")
                            result_items = self._default_value[name]

                result[name] = result_items
            # Dict
            elif att_type is dict :
                txt = widget.toPlainText().strip()
                d: Dict[str, Any] = {}
                if not txt == "" :
                    sample_v = next(iter(self._default_value[name].values()))
                    val_type = type(sample_v)
                    for line in txt.splitlines():
                        if ':' in line:
                            k, v = line.split(':', 1)
                            try:
                                d[k.strip()] = val_type(v.strip()) if v.strip() != '' else None
                            except Exception:
                                d[k.strip()] = v.strip()
                result[name] = d
            # fallback
            else:
                raise NotImplementedError(f"type {type(self._default_value[name])} not implemented for {name} parameters prompt.")
        return result