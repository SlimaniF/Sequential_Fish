"""
Submodule handling data opening and merging for viewer
"""
from types import NoneType, UnionType
from PyQt5.QtWidgets import (
    QPushButton, QHBoxLayout, QLabel,
    QDialog, QFormLayout, QLineEdit, QSpinBox, QDoubleSpinBox,QCheckBox, QWidget, QPlainTextEdit
)
from typing import Dict, Any,get_args, get_origin, Type
from pydantic import BaseModel


class ParametersModifier(QDialog):
    def __init__(self, data_model : Type[BaseModel], parent: QWidget | None = None, **parameters: Any):
        super().__init__(parent)
        self.setWindowTitle("Parameter Modification")
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
                default_value = {} if self._default_value[name] is None else self._default_value[name]
                lines = [f"{k}:{v}" for k, v in default_value.items()]
                pt.setPlainText("\n".join(lines)); pt.setFixedHeight(80)
                widget = pt

            # Fallback: show as string
            else:
                raise NotImplementedError(f"not gui implementation for type {name}; {att_type} --> {type_origin} from {default.annotation}. {att_args}")

            if add_none_option :
                sub_layout = QHBoxLayout()
                sub_layout.addWidget(widget)
                cb = QCheckBox()
                cb.setChecked(not self._default_value[name] is None)
                sub_layout.addWidget(cb)
                self._widgets[name+"__tuple__"] = [widget, cb]
                form.addRow(QLabel(name.replace('_', ' ').capitalize()), sub_layout)
    
            else :
                form.addRow(QLabel(name.replace('_', ' ').capitalize()), widget)
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

<<<<<<< HEAD
    
=======
>>>>>>> main
    def get_parameters(self) -> Dict[str,Any] :

        result: Dict[str, Any] = {}
        for name, default in self._param_types.items():
            att_type = default.annotation
            if not name in self._widgets.keys() or get_origin(att_type) is tuple:
                name += "__tuple__"

            widget = self._widgets[name]

            if get_origin(att_type) is UnionType:
                allowed_types = list(get_args(default.annotation))
                
                if NoneType in allowed_types :
                    allowed_types.remove(NoneType)
                else :
                    raise AssertionError("UnionType is allowed only to allow for NoneType.")
                assert len(allowed_types) == 1, allowed_types
                att_type = allowed_types[0]

                if not self._widgets[name + "__tuple__"][1].isChecked() : #User unselected parameter
                    result[name] = None
                    continue

            try :
                if "__tuple__" in name : name = name.replace("__tuple__","")
                result[name] = convert_input_to_type(widget, att_type)

            except InputError as e:
                raise InputError(f"Input error for : {name}") from e

            except Exception as e :
                raise ValueError(f"Critical error for {name}") from e

        return result


def convert_input_to_type(_input, to_type) :

    res = None
    allow_None = False
    if get_origin(to_type) is UnionType:

        allow_None = True
        allowed_types = list(get_args(to_type))
        
        if NoneType in allowed_types :
            allowed_types.remove(NoneType)
        else :
            raise AssertionError("UnionType is allowed only to allow for NoneType.")
        assert len(allowed_types) == 1, allowed_types
        to_type = allowed_types[0]


    if get_origin(to_type) is list :

        if not isinstance(_input, str) : _input = _input.text()
        if allow_None and _input == "" : return None
        
        allow_None_elmts = False
        type_args = list(get_args(to_type))

        if NoneType in type_args :
            allow_None_elmts = True
            if len(type_args) > 2 :raise ValueError(f"Expected 2 types in UnionType, a classic type and NoneType, got {type_args}")
            type_args.remove(NoneType)
            assert len(type_args) == 1, type_args
            elmts_type = type_args[0]
        else :
            assert len(type_args) == 1, type_args
            elmts_type = type_args[0]

        res = _process_list(_input, elmts_type=elmts_type, allow_None=allow_None_elmts)

    elif get_origin(to_type) is dict :
        
        if allow_None and _input == "" : return None

        if not isinstance(_input, str) : _input = _input.toPlainText().strip()
        allow_None_elmts = False

        key_type, values_type = get_args(to_type)
        if key_type not in [int, str] : raise TypeError(f"Supported types for dict keys are int or str not {key_type}")

        if get_origin(values_type) is UnionType :
            allowed_types = list(get_origin(values_type))
            if len(allowed_types) > 2 :raise ValueError(f"Expected 2 types in UnionType, a classic type and NoneType, got {allowed_types}")

            if NoneType in allowed_types :
                allowed_types.remove(NoneType)
                allow_None_elmts = True
                values_type = allowed_types[0]
            else :
                raise ValueError(f"Expected 2 types in UnionType, a classic type and NoneType, got no NoneType : {allowed_types}")

        res = _process_dict(_input, key_type, values_type, allow_None=allow_None_elmts)

    elif get_origin(to_type) is tuple :
        
        if isinstance(_input, str) : raise AssertionError("tuples types cannot be imbricated in another type such as list or dict.")
        
        elmt_types = get_args(to_type)

        if all(x is int for x in elmt_types):
            entries = tuple(sb.value() for sb in _input)
        elif all(x is str for x in elmt_types):
            entries = tuple(le.text() for le in _input)
        else : 
            raise NotImplementedError()
        
        res = _process_tuple(entries, elmts_types=elmt_types)
        
    elif to_type is int :
        if allow_None and _input == "" : return None

        if not isinstance(_input, (str,int)) : _input = _input.value()
        res = int(_input)
    elif to_type is str :
        if allow_None and _input == "" : return None
        if not isinstance(_input, (str,int,float)) : _input = _input.text()
        res = str(_input)
    elif to_type is float :
        if allow_None and _input == "" : return None
        if not isinstance(_input, (str,int,float)) : _input = _input.value()
        res = float(_input)
    elif to_type is bool :
        if isinstance(_input, str) :
            if _input in ["1","True","TRUE","true"] :
                res = True
            else :
                res = False
        else :
            res = bool(_input.isChecked())
    else :
        raise NotImplementedError(f"type {to_type} not implemented for {_input} parameters prompt.")
    
    return res

def _process_dict(dict_string : str, key_type, values_type, allow_None : bool) -> dict :

    res = {}
    for entry in dict_string.split('\n') :
        if ':' not in entry : raise InputError(f"Dict type was wrongly set by user, each entry expects a ':' to seperate key and value but no ':' could be find in {entry}")
        
        splitted_entry = entry.split(":")
        if len(splitted_entry) != 2 : raise InputError(f"Dict type was wrongly set by user, each entry expects a ':' to seperate key and value more than one ({len(splitted_entry)}) ':' was found.")
        key = key_type(splitted_entry[0].strip())
        val = splitted_entry[1].strip()
        if allow_None and val == "" :
            val = None
        else :
            val = convert_input_to_type(val, to_type=values_type)

        res[key] = val

    return res

def _process_list(list_string : str, elmts_type, allow_None : bool) -> list:

    if not isinstance(list_string, str) : raise TypeError(f"Expected string type to convert to list, got : {type(list_string)}")

    list_string = list_string.strip()

    if "," in list_string :
        res = []
        for elmt in list_string.split(',') :
            res.append(convert_input_to_type(elmt.strip(), elmts_type))
    else :
        res = [convert_input_to_type(list_string.strip(), elmts_type)]
    
    if allow_None :
        for i in range(len(res)) :
            if res[i] == "" : res[i] = None

    return res

def _process_tuple(list_elmts : tuple, elmts_types : tuple) -> tuple :
    
    if len(list_elmts) != len(list_elmts) :
        raise ValueError(f"Expected {len(elmts_types)} elmts from type declaration, found {len(list_elmts)} elements.")

    if NoneType in elmts_types : raise TypeError("NoneType is not allowed in tuples")
    if UnionType in elmts_types : raise TypeError("UnionType is not allowed in tuples")
    
    res = []
    for elmt, to_type in zip(list_elmts, elmts_types) :
        res.append(convert_input_to_type(elmt, to_type))
    
    res = tuple(res)
    return res

class InputError(Exception) :
    """
    Execption raised when user didn't respect input rules when settings parameters, this exception should not be caused a failure in by code logic.
    """